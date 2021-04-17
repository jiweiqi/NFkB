function getRateParams()



    ## initalize 
    n = zeros(69, 1); # NFkB Activation Module
    i = zeros(34, 1); # IKK Activation Module
    a = zeros(32, 1);
    t = zeros(17, 1);  # TNF production module 
    
    ## ----- IkB:NFkB Module [including A20 synthesis/deg] -----
    #           ikba       ikbb        ikbe
    n[1:3]  = [7e-5      1e-5       1e-6 ]'; # constitutive txn
    n[4:6]	= [8         0.02       0.3  ]'; # inducible txn
    n[7:9]  = [4.159060426         3          4.999962344    ]'; # Hill Coefficient
    n[10:12] = [0         37         37   ]'; # inducible txn delay
    n[13] = 0.005830947	; # mRNA Degradation [estimated in Lee et al. 2018 Processes]
    n[14:15] = [      3e-3       4e-3 ]'; # mRNA Degradation
    n[16:18] = [0.25      0.25       0.25 ]'; # translation rate
    
    n[19:21] = [0.09      9e-3       0.045]'; # Free IkB Import
    n[23:25] = [0.012     0.012      0.012]'; # Free IkB Export
    n[27:29] = [0.276     0.0276     0.138]'; # IkB:NFkB Import
    n[30:32] = [0.828     0.414      0.414]'; # IkB:NFkB Export
    n[22]   = 5.4;   # Free NFkB import()
    n[26]   = 0.0048;# Free NFkB export
    
    n[33:35] = [0.12      0.18       0.18 ]'; # IkB deg cytoplasm
    n[36:38] = [0.12      0.18       0.18 ]'; # IkB deg nucleus
    n[39:41] = [6e-5      6e-5       6e-5 ]'; # Bound IkB deg cyt
    n[42:44] = [6e-5      6e-5       6e-5 ]'; # Bound IkB deg nuc
    
    n[45:47] = [30        30         30   ]'; # IkB:NFkB Asn cyt
    n[48:50] = [30        30         30   ]'; # IkB:NFkB Asn nuc
    n[51:53] = [6e-5      6e-5       6e-5 ]'; # IkB:NFkB Dsn cyt
    n[54:56] = [6e-5      6e-5       6e-5 ]'; # IkB:NFkB Dsn nuc
    
    n[57:59] = [0.36      0.12       0.18 ]'; # Free IkB + IKK deg
    n[60:62] = [0.36      0.12       0.18 ]'; # IkB:NFkB + IKK deg
    
    # A20 [NEW]
    n[63]  = 2e-6;    # basal txn rate 2e-6
    n[64]  = 0.4;     # inducible txn for A[20],0.4
    n[65]  = 3;       # Hill Coefficient
    n[66]  = 0.035;   # mRNA Degradation
    n[67]  = 0.25;    # translation rate
    n[68]  = 0.0029;  # protein degradation
    n[69]  = 120;     # promoter shutdown [experiments show ~120min]
    
    ## ----- IKK Activation Module -----
    i[1] =     0.1698;  # k_b_lps 6
    i[2] =    0.17765;  # k_in_lps 7 
    i[3] =    0.26114 ;  # k_out_lps 8
    i[4] =     13.367;  # k_d_lpsen, 9
    
    # ligand binding
    i[5]   = 0.19 ;   # k_b_tlr4lps, 10
    i[6]   = 2.7  ;   # k_ub_tlr4lps, 11 
    i[7]   = 0.19 ;   # k_b_tlr4lpsen, 12
    i[8]   = 2.7  ;   # k_ub_tlr4lpsen, 13 
    
    # generation & degradation
    i[9]   =   0.037479238 ;   # k_g_tlr4, 14 estimated in Lee et al. 2018
    i[10]  =    0.89603;   # k_d_tlr4,  15
    i[11]  =    9036.88927	;   # k_d_tlr4en, 16 estimated in Lee et al. 2018
    i[16]  =     14.414;     # k_deg_TLR4LPS,  21
    i[17]  =    0.42132;    # k_d_tlr4lpsen, 22
    
    # shutting
    i[12]  =     0.1344  ;   # k_in_tlr4, 17
    i[13]  =     3.6099  ;   # k_out_tlf4, 18
    i[14]  =   0.23513  ;     # k_in_tlr4lps, 19
    i[15]  =   0.041496  ;   # k_out_tlr4lps, 20 
    
    # activation
    i[18]  =   3.2907 ;    # k_f_myd88 TL activated MyD88 -> (MyD88)_6, 23
    i[19]  = 3; # Hill coefficient for MyD88 activation, 24
    i[20]  =   0.057936;      #  km_a_myd88, 25
    i[21]  =    0.27504;     # k_dis_myd88, 0.02, 26
    i[22]  =    0.38509  ;     # k_a_trif  1e+8, 27
    i[23]  =   0.011822   ;     # k_i_trif, 0.01, 28
    
    # TRAF6
    i[24]  =     7.477 ; # k_a_TRAF6_MyS88s 66
    i[25]  =     3.4132; # k_a_TRAF6_TRIFs 67
    i[26]  =    0.21778;       # k_i_TRAF6  68 
                           
    # IKKK
    i[27]  =    0.34317;   # IKKK_off --> IKKK [TRAF6s mediated], 71
    i[28]  =    5e-7;   # IKKK_off --> IKKK [constitutive],72
    i[29]  = 0.25;    # IKKK     --> IKKK_off [constitutive], 73
    
    # IKK
    i[30]  = 4753.870211;    # IKK_off  --> IKK [IKKK mediated],74 [estimated in Lee et al. 2018]
    i[31]  = 5e-5 ;    # IKK_off  --> IKK [constitutive],75
    i[32]  = 0.02  ;    # IKK      --> IKK_off [constitutive], 76
    i[33]  = 0.028466263 ;    # IKK      --> IKK_i [constitutive]  , 77 [estimated in Lee et al. 2018]
    i[34]  = 0.02  ;    # IKKi     --> IKK_off [constitutive] , 78
    
    
    ## the TNFR parameters
    # ----- IKK Activation Module ----- NEW f
    a[1]   = 0.0154; # pd_m_tnf 45' half life of exogenous TNF, 49
    
    # tnfrm metabolism [synthesis & degradation]
    a[2]   = 2e-7;   # tnfrm synthesis [txn, tsl, localization], 36
    a[3]   = 0.0058; # tnfrm --> deg  -- 120' halflife, 37
    
    # TNF-Independent C1 Activation
    a[4]   = 1e-5;   # 3tnfrm --> TNFR, 38
    a[5]   = 0.1;    # TNFR   --> 3tnfrm, 39
    a[6]   = 0.023;  # TNFR internalization -- 30' halflife, 40 
    
    a[7]   = 100;    # TNFR + TTR --> C1_off, 41
    a[8]   = 0.1;    # C1_off --> TNFR + TTR, 42
    a[9]   = 30;     # C1_off --> C1, 43
    a[10]  = 2;      # C1     --> C1_off, 44
    a[11]  = 1000;   # C1     --> C1_off [A20 Mediated], 45
    a[12]  = 0.1;   # C1     --> TNFR + TTR, 46
    a[13]  = 0.023;   # C1_off internalization, 47
    a[14]  = 0.023;   # C1 internalization, 48
    
    # TNF-dependent C1 Activation
    a[15]  = 1100;   # 3tnfrm + tnf --> TNFRtnf, 50
    a[16]  = 1100;  # TNFR + tnf --> TNFRtnf, 51
    a[17]  = 0.021;  # TNFRtnf   --> TNFR + tnf, 52
    a[18]  = 0.023;   # TNFRtnf internalization -- 30' halflife, 53
    
    a[19]  = 100;   # TNFRtnf + TTR --> C1tnf_off, 54
    a[20]  = 0.1;   # C1tnf_off --> TNFRtnf + TTR, 55
    a[21]  = 30;   # C1tnf_off --> C1tnf, 56
    a[22]  = 2;  # C1tnf     --> C1tnf_off, 57
    a[23]  = 1000;  # C1tnf     --> C1tnf_off [A20 Mediated], 58
    a[24]  = 0.1;   # C1tnf     --> TNFRtnf + TTR, 59
    a[25]  = 0.023;   # C1tnf_off internalization, 60 
    a[26]  = 0.023;   # C1tnf internalization, 61
    
    a[27]  = 0.021;  # C1tnf_off --> C1_off + tnf, 62
    a[28]  = 1100;  # C1_off + tnf --> C1tnf_off, 63
    a[29]  = 0.021;  # C1tnf    --> C1 + tnf, 64
    a[30]  = 1100;  # C1 + tnf --> C1tnf, 65
    
    # IKKK
    
    a[31]  = 500;    # IKKK_off --> IKKK [C1 mediated]?500, 69
    a[32]  = 500;  # IKKK_off --> IKKK [C1tnf mediated], 70 
    
    

    ## TNF part
    t[1] = 0;# 1e-5;#tnf constitutive txn
    t[2] = 1e-5;# tnf induced txn [rescaled in the whole model], 150 
    t[3] = 2;# tnf transcription,Hill coefficient, 151
    t[4] = .65;# tnf transcription induction EC50, 152
    t[5] = .02;# tnf transcript deg. max rate, 154
    t[6] = .4;# tnf nascent process rate, 153
    t[7] = .05;# tnf protein syns rate, 156
    t[8] = .07;# tnf deg rate, 157
    t[9] = .07;# tnf secretion rate, 158
    
    # from private communication with Dr. Zhang Cheng in Prof. Hoffmann's group
    t[10] = 0.0001;  # K_a0
    t[11] = 0.1;     # K_a
    t[12] = 0.4;     # K_i
    
    
    t[13] = 833.3; # estimated
    t[14] = 90;  # assumed time needed for full effects of BFA
    t[15] = 0.01; # assumed rate of deactivation of TRAF6* due to A20
     
    ##
    t[16] = 0.998934823;  # Effect from BFA on translation coefficient
    t[17] = 0.2;  # denominator term in translation inhibition

    return n, i, a, t
end

NP, IP, AP, TP = getRateParams()

mutable struct RateParams
    NP
    IP
    AP
    TP
    IKK_TOTAL
    flag_noTnfFeedback
end

IKK_TOTAL = sum(INIT[33:35])
flag_noTnfFeedback = false

v = RateParams(NP, IP, AP, TP, IKK_TOTAL, flag_noTnfFeedback)