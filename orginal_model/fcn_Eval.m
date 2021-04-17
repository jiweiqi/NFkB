function fcn_Eval()
[v.NP,v.IP,v.AP,v.TP] = getRateParams();  % call reaction parameters
load('INIT.mat') % call initial conditions of states (stored as a variable 'INIT')

v.IKK_TOTAL  = sum(INIT(33:35));

%% experimental data with LPS alone
  IkBa_exp=[1	0.750059394	0.732550197	0.664902215	0.539958209	0.7147038	1.18817208	1.247760483;
              1	0.671638839	0.537731722	0.356330319	0.472835442	0.614109941	1.280911298	1.348835829;
              1	0.606845246	0.293407607	0.148004809	0.546927401	0.456998625	1.188020406	1.146866781];
  %% standard error in this dataset
  IkBa_e=[0	0.121567582	0.131557583	0.109118139	0.09695737	0.097321806	0.103929938	0.221212221;
      0	0.089239571	0.124826811	0.107842155	0.100927013	0.121127315	0.136760739	0.262838186;
      0	0.116947655	0.083573102	0.035121146	0.13506598	0.108440483	0.17673781	0.162245212];
  
%% LPS injection without BFA
    v.flag_noTnfFeedback=0; % indicator that BFA is not added
    dose=[ 0.001 0.005 0.025 ].*v.TP(13); % four different input conditions
    t_sam=[0 10 20 30 60 120 240 360]; % time instants when the measurements are taken

    options = odeset('RelTol', 1e-8,'AbsTol',1e-8,'NonNegative',1:length(INIT));
    titstr={'10 ng/mL','50 ng/mL','250 ng/mL'};
    for i=1:length(dose)
         INIT(18)=dose(i);

         [t,y]=ode15s(@(t,x) nfkbOde(t,x,v),[0:0.1:t_sam(end)],INIT,options);
         
         IkBat=sum(y(:,1:4),2);
         IkBat_fc(:,i)=(IkBat)/IkBat(1);  % first measurement
        Titstr=strcat('LPS concentration = ',titstr{i},' without BFA');
        figure(i)
        plot(t,IkBat_fc(:,i),'-')
        hold on
        errorbar(t_sam,IkBa_exp(i,:),IkBa_e(i,:),'o')
        hold off
        xlabel('Time (min)')
        ylabel('Fold Changes in I\kappa B\alpha Expression')
        title(Titstr)
         
    end
    num_fig=i;


%% experimental data with LPS+BFA
    TNF_exp=[1	1.50905227	1.151694217	0.997088604	1.111627033	1.144707201	2.420483244	2.658460625;
             1	0.918739139	0.966665332	0.92533487	1.222323374	10.0196476	42.58450406	59.56395661 ;
              1	0.798976592	0.88693086	0.897727998	3.264823692	23.72754934	56.35665097	70.03515026;
             1	1.082706458	0.738792823	0.895783662	8.031737518	30.76604615	51.14124965	58.85066254];
  
    IkB_GO_exp=[1	0.793306935	0.6750918	0.680469971	0.674445238	0.543568457	0.440558657	0.411863461;
                1	0.810821933	0.804665267	0.741895289	0.521451662	0.634390948	0.654101829	0.500794234;
                1	0.829076692	0.683837537	0.463205155	0.420849402	0.564274738	0.754542786	0.618588695;
               1	0.574545251	0.301445775	0.157461953	0.4126428	0.382258542	0.426676215	0.516500755];
  
 %% standard error in the experimental data (LPS+BFA)          
        TNF_e=[0	0.393556355	0.299819139	0.183874171	0.111023204	0.105685897	0.734633812	0.728543936;
               0	0.160944611	0.12032246	0.105079274	0.095555912	2.035086826	10.04437486	11.13787838;
               0	0.044966985	0.030100954	0.047649911	0.473839718	2.566887928	7.385683273	13.79880433;
               0	0.104499708	0.086911164	0.083377626	1.57604415	2.211378941	8.045121393	11.72785711];
           
        IkBa_e=[0	0.206420367	0.111336721	0.11856887	0.041249194	0.017103163	0.018819203	0.003299805;
                0	0.092979421	0.089878956	0.08753009	0.066390394	0.072594344	0.088740348	0.090023262;
               0	0.091200332	0.092643138	0.077278932	0.108585775	0.109907777	0.186857827	0.118242834;
               0	0.11764559	0.069675449	0.038525682	0.084284143	0.099460942	0.098673324	0.126125962];
           
    %% LPS injection along with BFA
    v.flag_noTnfFeedback=1; % indicator that BFA is added
    dose=[0 0.001 0.005 0.025 ].*v.TP(13); % four different input conditions
    t_sam=[0 10 20 30 60 120 240 360]; % time instants when the measurements are taken
    options = odeset('RelTol', 1e-8,'AbsTol',1e-8,'NonNegative',1:length(INIT));
    titstr={'0','10 ng/mL','50 ng/mL','250 ng/mL'};
    for i=1:length(dose)
         INIT(18)=dose(i);

         [t,y]=ode15s(@(t,x) nfkbOde(t,x,v),[0:0.1:t_sam(end)],INIT,options);
         
         IkBat=sum(y(:,1:4),2);
         IkBat_fc(:,i)=(IkBat)/IkBat(1);  % first measurement
         TNF_fold(:,i)=(y(:,38))/y(1,38);  % second measurement
        Titstr=strcat('LPS concentration = ',titstr{i},' along with BFA');
        figure(2*i-1+num_fig)
        plot(t,IkBat_fc(:,i),'-')
        hold on
        errorbar(t_sam,IkB_GO_exp(i,:),IkBa_e(i,:),'o')
        hold off
        xlabel('Time (min)')
        ylabel('Fold Changes in I\kappa B\alpha Expression')
        title(Titstr)
         figure(2*i+num_fig)
        plot(t,TNF_fold(:,i),'-')
        hold on
        errorbar(t_sam,TNF_exp(i,:),TNF_e(i,:),'o')
        hold off
        xlabel('Time (min)')
        ylabel('Fold Changes in TNF\alpha Expression')
        title(Titstr)
         
    end
   


end
