using Plots
using DifferentialEquations
using DelimitedFiles
using MAT
using ForwardDiff
using Zygote
using Flux
using ProgressBars, Printf
using Flux.Optimise:update!
using Flux.Losses:mae, mse

exp_data = matread("./orginal_model/exp_data.mat");
TNF_exp = exp_data["TNF_exp"];
TNF_e = exp_data["TNF_e"];
IkBa_exp = exp_data["IkBa_exp"];
IkBa_e = exp_data["IkBa_e"];
IkB_GO_exp = exp_data["IkB_GO_exp"];
IkB_GO_e = exp_data["IkB_GO_e"];
INIT = exp_data["INIT"];

include("getRateParams.jl")
include("nfkbOde_clean.jl")

## LPS injection without BFA
v.flag_noTnfFeedback = false; # indicator that BFA is not added
dose = vec([ 0.001 0.005 0.025 ] .* v.TP[13]); # four different input conditions
t_sam = vec(Float64.([0 10 20 30 60 120 240 360])); # time instants when the measurements are taken

## Test AD
i = 1;
INIT[18] = dose[i]
u0 = vec(INIT);
delta = similar(u0);
t = 0.01;
nfkbOde!(delta, u0, v, t)
gu = ForwardDiff.jacobian((du, u) -> nfkbOde!(du, u, v, t), delta, u0)

n_nodes = 2;
ns = length(u0);
nn = Chain(x -> x,
        Dense(ns, n_nodes, gelu),
        Dense(n_nodes, n_nodes, gelu),
        Dense(n_nodes, n_nodes, gelu),
        Dense(n_nodes, n_nodes, gelu),
        Dense(n_nodes, ns))
p, re = Flux.destructure(nn);

function hybrid!(delta, u, p, t)
    nfkbOde!(delta, u, v, t);
    du = re(p)(u)
    delta .+= du
end
hybrid!(delta, u0, p, t)


titstr = ["10, ng/mL", "50, ng/mL", "250, ng/mL"];

for i = 1:length(dose)
    INIT[18] = dose[i]

    u0 = vec(INIT);
    prob = ODEProblem(hybrid!, u0, (0.0, t_sam[end]), p);
    sol = solve(prob, atol=1.e-6, rtol=1.e-3);
    IkBat = sum(sol'[:,1:4], dims=2)
    
    plt = plot(sol.t, IkBat ./ IkBat[1], lw=3);
    scatter!(plt, t_sam, IkBa_exp[i,:], yerr=IkBa_e[i,:]);
    xlabel!("Time [min]");
    ylabel!("Fold Changes in Iκ Bα Expression");
    title!(string("LPS concentration = ", "$(titstr[i])", " without BFA"));
    png(plt, "./figs/IkB_$i without BFA.png");
end

## LPS injection along with BFA
v.flag_noTnfFeedback = true; # indicator that BFA is added
titstr = ["0", "10, ng/mL", "50, ng/mL", "250, ng/mL"];

for i = 1:length(dose)
    INIT[18] = dose[i]

    u0 = vec(INIT);
    prob = ODEProblem(nfkbOde!, u0, (0.0, t_sam[end]), v);
    sol = solve(prob, atol=1.e-6, rtol=1.e-3);
    IkBat = sum(sol'[:,1:4], dims=2)
    
    plt = plot(sol.t, IkBat ./ IkBat[1], lw=3);
    scatter!(plt, t_sam, IkB_GO_exp[i,:], yerr=IkB_GO_e[i,:]);
    xlabel!("Time [min]");
    ylabel!("Fold Changes in Iκ Bα Expression");
    title!(string("LPS concentration = ", "$(titstr[i])", " with BFA"));
    png(plt, "./figs/IkB_$i with BFA.png");

    plt = plot(sol.t, sol'[:, 38] ./ sol'[1, 38], lw=3);
    scatter!(plt, t_sam, TNF_exp[i,:], yerr=TNF_e[i,:]);
    xlabel!("Time [min]");
    ylabel!("Fold Changes in TNFα Expression");
    title!(string("LPS concentration = ", "$(titstr[i])", " with BFA"));
    png(plt, "./figs/TNF_$i with BFA.png");
end