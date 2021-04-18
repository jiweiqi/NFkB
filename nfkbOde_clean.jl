function nfkbOde!(delta, x, v, t)

    # nf-kb part
    IkBa           = x[1]
    IkBan          = x[2]
    IkBaNFkB       = x[3]
    IkBaNFkBn      = x[4]
    IkBat          = x[5]
    IkBb           = x[6]
    IkBbn          = x[7]
    IkBbNFkB       = x[8]
    IkBbNFkBn      = x[9]
    IkBbt          = x[10]
    IkBe           = x[11]
    IkBen          = x[12]
    IkBeNFkB       = x[13]
    IkBeNFkBn      = x[14]
    IkBet          = x[15]
    NFkB           = x[16]
    NFkBn          = x[17]
    
    
    # TLR4&IKK part
    LPSo           = x[18]; # input LPS
    LPS            = x[19]; 
    LPSen          = x[20]; 
    TLR4           = x[21]; 
    TLR4en         = x[22];      
    TLR4LPS        = x[23]
    TLR4LPSen      = x[24];        
    MyD88          = x[25]
    MyD88s         = x[26]
    TRIF           = x[27]; 
    TRIFs          = x[28]; 
    TRAF6          = x[29]
    TRAF6s         = x[30]
    IKKK_off       = x[31]; 
    IKKK           = x[32]; 
    IKK_off        = x[33]; 
    IKK            = x[34]; 
    IKK_i          = x[35]; 
    
    
    # TNF part
    TNFnas        = x[36]
    TNFmRNA       = x[37]
    TNFpro        = x[38]
    TNF           = x[39]
    
    # TNFR part 
    tnfrm          = x[40]
    TNFR           = x[41]
    TNFRtnf        = x[42]
    C1             = x[43]
    C1_off         = x[44]
    C1tnf          = x[45]
    C1tnf_off      = x[46]
    TTR            = x[47]
    # A20 
    a20            = x[48]
    a20t           = x[49]

    ## set IKK flux
    IKK_flux = IKK / (v.IKK_TOTAL); # ikk modification
    
    # TNF mRNA stabilization 
    fa = (TRIFs + v.TP[10]) / (TRIFs + v.TP[11])
    fi = v.TP[12] / (TRIFs + v.TP[12])
    
    kdegm = v.TP[5] * fi
    kpr = v.TP[6] * fa
    ktl = v.TP[7] * fa
    ksec = v.TP[9] * fa
    
    # golgi plug effect
    kdeg_tnf = v.TP[8]
    ktl_a = v.NP[16]
    ktl_b = v.NP[17]
    ktl_e = v.NP[18]
     
    # TODO:  is t > 0 nessessary?
    if (v.flag_noTnfFeedback) & (t > 0.0)
        golgi = t / (t + v.TP[14])
        ksec = ksec * (1 - golgi)
        ktnfr = v.AP[2] * (1 - golgi)
        ktlr4 = v.IP[9] * (1 - golgi)
        kshut = v.IP[2] * (1 - golgi); 
        kshut_com = v.IP[14] * (1 - golgi)
    else
        ktnfr = v.AP[2]
        ktlr4 = v.IP[9]
        kshut = v.IP[2]; 
        kshut_com = v.IP[14]
    end
    
    ## Calculate Iteration Fluxes
    
    # -- NFkB Activation Module
    f_rsu_a         = v.NP[1] 
    f_rsu_b         = v.NP[2] 
    f_rsu_e         = v.NP[3] 
    f_rsu_f         = v.TP[1] ; # tnf basal transcription [NEW]
    
    f_rsr_an        = v.NP[4]   * (NFkBn^v.NP[7]) 
    f_rsr_bn        = v.NP[5]   * (NFkBn^v.NP[8]) 
    f_rsr_en        = v.NP[6]   * (NFkBn^v.NP[9]) 
    
    f_rsr_fn        = v.TP[2]  * (NFkBn^v.TP[3]) / (NFkBn^v.TP[3] + v.TP[4]^v.TP[3]); # tnf transcription 
    f_rp_f          = kpr * TNFnas; # tnf RNA processing
    
    f_rd_a          = v.NP[13]  * IkBat	
    f_rd_b          = v.NP[14]  * IkBbt	
    f_rd_e          = v.NP[15]  * IkBet	
    f_rd_f          = kdegm * TNFmRNA; # tnf mRNA degradation
    
    f_ps_c_a        = ktl_a	* IkBat	
    f_ps_c_b        = ktl_b	* IkBbt	
    f_ps_c_e        = ktl_e	* IkBet	
    f_ps_c_f        = ktl	* TNFmRNA;  # tnf protein synthsis 
    f_pd_f          = kdeg_tnf    * TNFpro;  # tnf protein deg. (inside the cell)
    f_sc_c_f        = ksec   * TNFpro;  # tnf protein secrection 
    
    f_in_a          = v.NP[19]	* IkBa	
    f_in_b          = v.NP[20]	* IkBb	
    f_in_e          = v.NP[21]	* IkBe	
    f_in_n          = v.NP[22]	* NFkB	
    f_ex_a          = v.NP[23]	* IkBan	
    f_ex_b          = v.NP[24]	* IkBbn	
    f_ex_e          = v.NP[25]	* IkBen	
    f_ex_n          = v.NP[26]	* NFkBn	
    f_in_2an        = v.NP[27]	* IkBaNFkB	
    f_in_2bn        = v.NP[28]	* IkBbNFkB	
    f_in_2en        = v.NP[29]	* IkBeNFkB	
    f_ex_2an        = v.NP[30]	* IkBaNFkBn	
    f_ex_2bn        = v.NP[31]	* IkBbNFkBn	
    f_ex_2en        = v.NP[32]	* IkBeNFkBn	
    
    f_pd_c_a        = v.NP[33]	* IkBa	
    f_pd_c_b        = v.NP[34]	* IkBb	
    f_pd_c_e        = v.NP[35]	* IkBe	
    f_pd_n_a        = v.NP[36]	* IkBan	
    f_pd_n_b        = v.NP[37]	* IkBbn	
    f_pd_n_e        = v.NP[38]	* IkBen	
    f_pd_c_2an      = v.NP[39]	* IkBaNFkB	
    f_pd_c_2bn      = v.NP[40]	* IkBbNFkB	
    f_pd_c_2en      = v.NP[41]	* IkBeNFkB	
    f_pd_n_2an      = v.NP[42]	* IkBaNFkBn	
    f_pd_n_2bn      = v.NP[43]	* IkBbNFkBn	
    f_pd_n_2en      = v.NP[44]	* IkBeNFkBn	
    
    f_a_c_an        = v.NP[45]  * IkBa	* NFkB	
    f_a_c_bn        = v.NP[46]  * IkBb	* NFkB	
    f_a_c_en        = v.NP[47]  * IkBe	* NFkB	
    f_a_n_an        = v.NP[48]	* IkBan	* NFkBn	
    f_a_n_bn        = v.NP[49]	* IkBbn	* NFkBn	
    f_a_n_en        = v.NP[50]  * IkBen	* NFkBn	
    f_d_c_an        = v.NP[51]	* IkBaNFkB	
    f_d_c_bn        = v.NP[52]	* IkBbNFkB	
    f_d_c_en        = v.NP[53]	* IkBeNFkB	
    f_d_n_an        = v.NP[54]	* IkBaNFkBn	
    f_d_n_bn        = v.NP[55]	* IkBbNFkBn	
    f_d_n_en        = v.NP[56]	* IkBeNFkBn	
    
    # IKK Mediated IkB Degradation [free & bound]
    f_ph_c_a        = v.NP[57]  * IkBa * IKK_flux 
    f_ph_c_b        = v.NP[58]	* IkBb * IKK_flux 
    f_ph_c_e        = v.NP[59]	* IkBe * IKK_flux 
    
    f_ph_c_an       = v.NP[60]  * IkBaNFkB * IKK_flux 
    f_ph_c_bn       = v.NP[61]	* IkBbNFkB * IKK_flux 
    f_ph_c_en       = v.NP[62]	* IkBeNFkB * IKK_flux 
    
    # A20 Fluxes 
    f_rsu_a20       = v.NP[63] 
    f_rsr_a20       = v.NP[64]  * (NFkBn^v.NP[65]) 
    f_rd_a20        = v.NP[66]  * a20t 
    f_ps_c_a20      = v.NP[67]  * NFkBn 
    f_pd_c_a20      = v.NP[68]  * a20  
    
    if t > v.NP[69]  # disables inducible A20 txn to match exp data
        f_rsr_a20 *= 0.0
    end
    
    # -- IKK Activation Module & TLR4 module
    # LPS binding to the surface()
    f_b_LPS         = v.IP[1]   * LPSo  
    
    # ligand binding
    f_b_TLR4LPS     = v.IP[5]    * TLR4 * LPS 
    f_ub_TLR4LPS    = v.IP[6]    * TLR4LPS 
    f_b_TLR4LPSen   = v.IP[7]    * TLR4en * LPSen 
    f_ub_TLR4LPSen  = v.IP[8]    * TLR4LPSen 
    
    # Activation modified on 0426
    f_a_MyD88       = v.IP[18]    * (TLR4LPS)^v.IP[19] / ((TLR4LPS)^v.IP[19] + (v.IP[20])^v.IP[19])  * MyD88
    f_i_MyD88       = v.IP[21]    * MyD88s 
    f_a_TRIF        = v.IP[22]    * TLR4LPSen * TRIF
    f_i_TRIF        = v.IP[23]   * TRIFs 
    
    # shuttling; modified on 0426
    f_in_LPS        = kshut   * LPS 
    f_out_LPS       = v.IP[3]   * LPSen
    f_in_TLR4       = v.IP[12]   * TLR4
    f_out_TLR4      = v.IP[13]   * TLR4en
    f_in_TLR4LPS    = kshut_com   * TLR4LPS
    f_out_TLR4LPS   = v.IP[15]   * TLR4LPSen
    
    # generation & degradation
    f_g_TLR4        = ktlr4
    f_d_TLR4        = v.IP[10]   * TLR4
    f_d_TLR4en      = v.IP[11]   * TLR4en
    f_d_LPSen       = v.IP[4]   * LPSen
    f_d_TLR4LPSen   = v.IP[17]   * TLR4LPSen
    # modified 0515
    f_d_TLR4LPS     = v.IP[16]   * TLR4LPS
    # -----------------------IKK module
    # TRAF6
    f_a_TRAF6_MyD88s = v.IP[24]   * TRAF6 * MyD88s
    f_a_TRAF6_TRIFs = v.IP[25]   * TRAF6 * TRIFs
    f_i_TRAF6       = v.IP[26]   * TRAF6s
    
    # IKKK
    
    f_IKKK_on       = v.IP[28]   * IKKK_off
    f_IKKK_on_TRAF6s = v.IP[27]   * IKKK_off * TRAF6s
    f_IKKK_off      = v.IP[29]   * IKKK
    
    # IKK
    f_IKK_on        = v.IP[31]   * IKK_off
    f_IKK_on_IKKK   = v.IP[30]   * IKK_off * IKKK
    f_IKK_off       = v.IP[32]   * IKK
    f_IKK_off_i     = v.IP[33]   * IKK
    f_IKKi          = v.IP[34]   * IKK_i
    
    ## TNFR part 
    f_pd_m_tnf      = v.AP[1]    * TNF ; # pd_m_tnf 45' half life of exogenous TNF
    f_syn_tnfrm     = ktnfr;  # tnfrm synthesis [txn, tsl, localization]
    f_pd_tnfrm      = v.AP[3]    * tnfrm;# tnfrm --> deg  -- 120' halflife
    f_a_tnfrm       = v.AP[4]    * tnfrm;# 3tnfrm --> TNFR
    f_d_TNFR        = v.AP[5]    * TNFR;# TNFR   --> 3tnfrm
    f_i_TNFR        = v.AP[6]    * TNFR;# TNFR internalization -- 30' halflife
    f_a_C1_off      = v.AP[7]    * TNFR * TTR;# TNFR + TTR --> C1_off
    f_d_C1_off      = v.AP[8]    * C1_off;# C1_off --> TNFR + TTR
    f_a_C1          = v.AP[9]    * C1_off;# C1_off --> C1
    f_C1_off        = v.AP[10]   * C1;# C1     --> C1_off
    f_C1_A20        = v.AP[11]   * C1 * a20;# C1     --> C1_off [A20 Mediated]
    f_d_C1          = v.AP[12]   * C1; # UPDATE # C1     --> TNFR + TTR
    f_i_C1_off      = v.AP[13]   * C1_off; # C1_off internalization
    f_i_C1          = v.AP[14]   * C1; # C1 internalization
    f_a_tnfrmtnf    = v.AP[15]   * tnfrm * TNF; # 3tnfrm + tnf --> TNFRtnf
    f_a_TNFRtnf     = v.AP[16]   * TNFR * TNF; # TNFR + tnf --> TNFRtnf
    f_d_TNFRtnf     = v.AP[17]   * TNFRtnf; # TNFRtnf   --> TNFR + tnf
    f_i_TNFRtnf     = v.AP[18]   * TNFRtnf;  # TNFRtnf internalization -- 30' halflife
    
    f_a_C1tnf_off   = v.AP[19]   * TNFRtnf * TTR; # TNFRtnf + TTR --> C1tnf_off
    f_d_C1tnf_off   = v.AP[20]   * C1tnf_off; # C1tnf_off --> TNFRtnf + TTR
    f_a_C1tnf       = v.AP[21]   * C1tnf_off; # C1tnf_off --> C1tnf
    f_C1tnf_off     = v.AP[22]   * C1tnf; # C1tnf     --> C1tnf_off
    f_C1tnf_A20     = v.AP[23]   * C1tnf * a20;  # C1tnf     --> C1tnf_off [A20 Mediated]
    f_d_C1tnf       = v.AP[24]   * C1tnf; # updated  # C1tnf     --> TNFRtnf + TTR
    f_i_C1tnf       = v.AP[25]   * C1tnf; # C1tnf internalization
    f_i_C1tnf_off   = v.AP[26]   * C1tnf_off;  # C1tnf_off internalization
    f_d_tnf_C1_off  = v.AP[27]   * C1tnf_off; # C1tnf_off --> C1_off + tnf
    f_a_tnf_C1_off  = v.AP[28]   * C1_off * TNF; # C1_off + tnf --> C1tnf_off
    f_d_tnf_C1      = v.AP[29]   * C1tnf; # C1tnf    --> C1 + tnf
    f_a_tnf_C1      = v.AP[30]   * C1 * TNF; # C1 + tnf --> C1tnf
    f_IKKK_on_C1    = v.AP[31]   * IKKK_off * C1; # IKKK_off --> IKKK [C1 mediated]?500
    f_IKKK_on_C1tnf = v.AP[32]   * IKKK_off * C1tnf;# IKKK_off --> IKKK [C1tnf mediated]
    
    ## Add Fluxes to component concentrations; then Save
    
    # IkB Transcrv.IPtion
    d_IkBat     = d_IkBat + f_rsu_a
    d_IkBbt     = d_IkBbt + f_rsu_b
    d_IkBet     = d_IkBet + f_rsu_e
    d_IkBat     = d_IkBat + f_rsr_an
    d_IkBbt     = d_IkBbt + f_rsr_bn
    d_IkBet     = d_IkBet + f_rsr_en
    d_TNFnas    = d_TNFnas + f_rsu_f; # NEW
    d_TNFnas    = d_TNFnas + f_rsr_fn; # NEW
    d_TNFnas    = d_TNFnas - f_rp_f; # NEW
    d_TNFmRNA   = d_TNFmRNA + f_rp_f; # NEW
    
    
    # IkB Transcrv.IPt Degradation
    d_IkBbt     = d_IkBbt - f_rd_b
    d_IkBet     = d_IkBet - f_rd_e
    d_IkBat     = d_IkBat - f_rd_a
    d_TNFmRNA   = d_TNFmRNA  - f_rd_f;   # tnft degradation
    
    # IkB Translation
    d_IkBat     = d_IkBat - f_ps_c_a
    d_IkBat     = d_IkBat + f_ps_c_a
    d_IkBa      = + f_ps_c_a
    
    d_IkBbt     = d_IkBbt - f_ps_c_b
    d_IkBbt     = d_IkBbt + f_ps_c_b
    d_IkBb      = d_IkBb  + f_ps_c_b
    
    d_IkBet     = d_IkBet - f_ps_c_e
    d_IkBet     = d_IkBet + f_ps_c_e
    d_IkBe      = d_IkBe  + f_ps_c_e
    
    d_TNFpro    = d_TNFpro + f_ps_c_f; 
    
    # TNF secretion + protein degradation
    d_TNFpro     = d_TNFpro - f_sc_c_f; 
    d_TNF        = d_TNF    + f_sc_c_f; 
    d_TNFpro     = d_TNFpro - f_pd_f; 
    
    # IkB:NFkB Shuttling [Free & Bound]
    d_IkBa      = d_IkBa  - f_in_a
    d_IkBan     = + f_in_a
    
    d_IkBb      = d_IkBb  - f_in_b
    d_IkBbn     = d_IkBbn + f_in_b
    
    d_IkBe      = d_IkBe  - f_in_e
    d_IkBen     = d_IkBen + f_in_e
    
    d_NFkB      = d_NFkB  - f_in_n
    d_NFkBn     = d_NFkBn + f_in_n
    
    d_IkBan     = d_IkBan - f_ex_a
    d_IkBa      = d_IkBa  + f_ex_a
    
    d_IkBbn     = d_IkBbn - f_ex_b
    d_IkBb      = d_IkBb  + f_ex_b
    
    d_IkBen     = d_IkBen - f_ex_e
    d_IkBe      = d_IkBe  + f_ex_e
    
    d_NFkBn     = d_NFkBn - f_ex_n
    d_NFkB      = d_NFkB  + f_ex_n
    
    d_IkBaNFkB  = - f_in_2an
    d_IkBaNFkBn = d_IkBaNFkBn + f_in_2an
    
    d_IkBbNFkB  = d_IkBbNFkB  - f_in_2bn
    d_IkBbNFkBn = d_IkBbNFkBn + f_in_2bn
    
    d_IkBeNFkB  = d_IkBeNFkB  - f_in_2en
    d_IkBeNFkBn = d_IkBeNFkBn + f_in_2en
    
    d_IkBaNFkBn = d_IkBaNFkBn - f_ex_2an
    d_IkBaNFkB  = d_IkBaNFkB  + f_ex_2an
    
    d_IkBbNFkBn = d_IkBbNFkBn - f_ex_2bn
    d_IkBbNFkB  = d_IkBbNFkB  + f_ex_2bn
    
    d_IkBeNFkBn = d_IkBeNFkBn - f_ex_2en
    d_IkBeNFkB  = d_IkBeNFkB  + f_ex_2en
    
    # IkB:NFkB Association [Cytoplasm & Nucleus]
    d_IkBa      = d_IkBa - f_a_c_an
    d_NFkB      = d_NFkB - f_a_c_an
    d_IkBaNFkB  = d_IkBaNFkB + f_a_c_an
    
    d_IkBb      = d_IkBb - f_a_c_bn
    d_NFkB      = d_NFkB - f_a_c_bn
    d_IkBbNFkB  = d_IkBbNFkB + f_a_c_bn
    
    d_IkBe      = d_IkBe - f_a_c_en
    d_NFkB      = d_NFkB - f_a_c_en
    d_IkBeNFkB  = d_IkBeNFkB + f_a_c_en
    
    d_IkBan     = d_IkBan - f_a_n_an
    d_NFkBn     = d_NFkBn - f_a_n_an
    d_IkBaNFkBn = d_IkBaNFkBn + f_a_n_an
    
    d_IkBbn     = d_IkBbn - f_a_n_bn
    d_NFkBn     = d_NFkBn - f_a_n_bn
    d_IkBbNFkBn = d_IkBbNFkBn + f_a_n_bn
    
    d_IkBen     = d_IkBen - f_a_n_en
    d_NFkBn     = d_NFkBn - f_a_n_en
    d_IkBeNFkBn = d_IkBeNFkBn + f_a_n_en
    
    # IkB:NFkB Dissociation [Cytoplasm & Nucleus]
    d_IkBaNFkB  = d_IkBaNFkB - f_d_c_an
    d_IkBa      = d_IkBa + f_d_c_an
    d_NFkB      = d_NFkB + f_d_c_an
    
    d_IkBbNFkB  = d_IkBbNFkB - f_d_c_bn
    d_IkBb      = d_IkBb + f_d_c_bn
    d_NFkB      = d_NFkB + f_d_c_bn
    
    d_IkBeNFkB  = d_IkBeNFkB - f_d_c_en
    d_IkBe      = d_IkBe + f_d_c_en
    d_NFkB      = d_NFkB + f_d_c_en
    
    d_IkBaNFkBn = d_IkBaNFkBn - f_d_n_an
    d_IkBan     = d_IkBan + f_d_n_an
    d_NFkBn     = d_NFkBn + f_d_n_an
    
    d_IkBbNFkBn = d_IkBbNFkBn - f_d_n_bn
    d_IkBbn     = d_IkBbn + f_d_n_bn
    d_NFkBn     = d_NFkBn + f_d_n_bn
    
    d_IkBeNFkBn = d_IkBeNFkBn - f_d_n_en
    d_IkBen     = d_IkBen + f_d_n_en
    d_NFkBn     = d_NFkBn + f_d_n_en
    
    # Free IkB Degradation [Cytoplasm & Nucleus]
    d_IkBa      = d_IkBa  - f_pd_c_a
    d_IkBb      = d_IkBb  - f_pd_c_b
    d_IkBe      = d_IkBe  - f_pd_c_e
    d_IkBan     = d_IkBan - f_pd_n_a
    d_IkBbn     = d_IkBbn - f_pd_n_b
    d_IkBen     = d_IkBen - f_pd_n_e
    
    # IkB:NFkB Degradation [Cytoplasm & Nucleus]
    d_IkBaNFkB  = d_IkBaNFkB - f_pd_c_2an
    d_NFkB      = d_NFkB + f_pd_c_2an
    
    d_IkBbNFkB  = d_IkBbNFkB - f_pd_c_2bn
    d_NFkB      = d_NFkB + f_pd_c_2bn
    
    d_IkBeNFkB  = d_IkBeNFkB - f_pd_c_2en
    d_NFkB      = d_NFkB + f_pd_c_2en
    
    d_IkBaNFkBn = d_IkBaNFkBn - f_pd_n_2an
    d_NFkBn     = d_NFkBn + f_pd_n_2an
    
    d_IkBbNFkBn = d_IkBbNFkBn - f_pd_n_2bn
    d_NFkBn     = d_NFkBn + f_pd_n_2bn
    
    d_IkBeNFkBn = d_IkBeNFkBn - f_pd_n_2en
    d_NFkBn     = d_NFkBn + f_pd_n_2en
    
    # IKK Mediated IkB Degradation
    d_IkBa      = d_IkBa - f_ph_c_a
    d_IkBb      = d_IkBb - f_ph_c_b
    d_IkBe      = d_IkBe - f_ph_c_e
    
    # IKK Mediated IkB:NFkB Degradation
    d_IkBaNFkB  = d_IkBaNFkB - f_ph_c_an
    d_NFkB      = d_NFkB + f_ph_c_an
    
    d_IkBbNFkB  = d_IkBbNFkB - f_ph_c_bn
    d_NFkB      = d_NFkB + f_ph_c_bn
    
    d_IkBeNFkB  = d_IkBeNFkB - f_ph_c_en
    d_NFkB      = d_NFkB + f_ph_c_en
    
    
    # Upstream Pathway  [LPS] --> IKK activation
    
    # Ligand binding
    d_LPSo      = d_LPSo - f_b_LPS
    d_LPS       = d_LPS  + f_b_LPS
    
    d_LPS       = d_LPS       -   f_b_TLR4LPS
    d_TLR4       = d_TLR4       -   f_b_TLR4LPS
    d_TLR4LPS   = d_TLR4LPS   +   f_b_TLR4LPS
    d_LPS       = d_LPS       +   f_ub_TLR4LPS
    d_TLR4       = d_TLR4     +   f_ub_TLR4LPS
    d_TLR4LPS   = d_TLR4LPS   -   f_ub_TLR4LPS
    
    d_LPSen     = d_LPSen     -   f_b_TLR4LPSen
    d_TLR4en     = d_TLR4en     -   f_b_TLR4LPSen
    d_TLR4LPSen = d_TLR4LPSen +   f_b_TLR4LPSen
    d_LPSen     = d_LPSen     +   f_ub_TLR4LPSen
    d_TLR4en     = d_TLR4en   +   f_ub_TLR4LPSen
    d_TLR4LPSen = d_TLR4LPSen -   f_ub_TLR4LPSen
    
    # Activation
    d_MyD88     = d_MyD88     -   f_a_MyD88
    d_MyD88s    = d_MyD88s    +   f_a_MyD88
    d_MyD88     = d_MyD88     +   f_i_MyD88
    d_MyD88s    = d_MyD88s    -   f_i_MyD88
    
    d_TRIF      = d_TRIF      -   f_a_TRIF
    d_TRIFs     = d_TRIFs     +   f_a_TRIF
    d_TRIF      = d_TRIF      +   f_i_TRIF
    d_TRIFs     = d_TRIFs     -   f_i_TRIF
    
    # shuttling
    d_LPS       = d_LPS       -   f_in_LPS
    d_LPSen     = d_LPSen     +   f_in_LPS
    # d_LPSo       = d_LPSo     +   f_out_LPS
    d_LPS       = d_LPS       +   f_out_LPS
    d_LPSen     = d_LPSen     -   f_out_LPS
    
    d_TLR4      = d_TLR4      -   f_in_TLR4
    d_TLR4en    = d_TLR4en    +   f_in_TLR4
    d_TLR4      = d_TLR4      +   f_out_TLR4
    d_TLR4en    = d_TLR4en    -   f_out_TLR4
    
    d_TLR4LPS      = d_TLR4LPS      -   f_in_TLR4LPS
    d_TLR4LPSen    = d_TLR4LPSen    +   f_in_TLR4LPS
    d_TLR4LPS      = d_TLR4LPS      +   f_out_TLR4LPS
    d_TLR4LPSen     = d_TLR4LPSen   -   f_out_TLR4LPS
    
    
    # generation & degradation
    d_TLR4      = d_TLR4 + f_g_TLR4
    d_TLR4      = d_TLR4 - f_d_TLR4
    d_TLR4en     = d_TLR4en - f_d_TLR4en
    d_LPSen      = d_LPSen - f_d_LPSen
    d_TLR4LPSen      = d_TLR4LPSen - f_d_TLR4LPSen
    # -- modified 0515
    d_TLR4LPS      = d_TLR4LPS - f_d_TLR4LPS
    
    
    # TRAF6
    d_TRAF6       = d_TRAF6  - f_a_TRAF6_MyD88s
    d_TRAF6s      = d_TRAF6s + f_a_TRAF6_MyD88s
    d_TRAF6       = d_TRAF6  - f_a_TRAF6_TRIFs
    d_TRAF6s      = d_TRAF6s + f_a_TRAF6_TRIFs
    d_TRAF6       = d_TRAF6  + f_i_TRAF6
    d_TRAF6s      = d_TRAF6s - f_i_TRAF6
    
    # IKKK Regulation
    d_IKKK      = d_IKKK     + f_IKKK_on
    d_IKKK_off  = d_IKKK_off - f_IKKK_on
    
    d_IKKK      = d_IKKK     + f_IKKK_on_TRAF6s
    d_IKKK_off  = d_IKKK_off - f_IKKK_on_TRAF6s
    
    d_IKKK_off  = d_IKKK_off + f_IKKK_off
    d_IKKK      = d_IKKK     - f_IKKK_off
    
    # IKK Regulation
    d_IKK       = d_IKK     + f_IKK_on
    d_IKK_off   = d_IKK_off - f_IKK_on
    
    d_IKK       = d_IKK     + f_IKK_on_IKKK
    d_IKK_off   = d_IKK_off - f_IKK_on_IKKK
    
    d_IKK       = d_IKK     - f_IKK_off
    d_IKK_off   = d_IKK_off + f_IKK_off
    
    d_IKK       = d_IKK     - f_IKK_off_i
    d_IKK_i     = d_IKK_i   + f_IKK_off_i
    
    d_IKK_off   = d_IKK_off + f_IKKi
    d_IKK_i     = d_IKK_i   - f_IKKi
    
    ## TNR module [NEW] + A20 [NEW]
    # A20 Expression Transcrviption
    d_a20t      = d_a20t + f_rsu_a20
    d_a20t      = d_a20t + f_rsr_a20
    d_a20t      = d_a20t - f_rd_a20
    
    # Translation
    d_a20       = d_a20  + f_ps_c_a20
    
    # Degradation
    d_a20       = d_a20 - f_pd_c_a20
    
    # Upstream Pathway  [TNF] --> IKK activation
    
    # TNF-independent Activation
    
    # tnfrm trimerization 3tnfrm<--> TNFR
    d_tnfrm     = d_tnfrm - 3 * f_a_tnfrm
    d_TNFR      = d_TNFR +   f_a_tnfrm
    
    d_tnfrm     = d_tnfrm + 3 * f_d_TNFR
    d_TNFR      = d_TNFR -   f_d_TNFR
    
    # TNF Receptor Metabolism
    d_TNFR     = d_TNFR     - f_i_TNFR
    d_TNFRtnf  = d_TNFRtnf  - f_i_TNFRtnf
    
    d_tnfrm    = d_tnfrm    - f_pd_tnfrm
    d_tnfrm    = d_tnfrm    + f_syn_tnfrm
    
    # Complex I Generation  TTR + TNFR <--> C1_off; C1
    d_C1_off    = d_C1_off  + f_a_C1_off
    d_TTR       = d_TTR     - f_a_C1_off
    d_TNFR      = d_TNFR    - f_a_C1_off
    
    d_C1_off    = d_C1_off  - f_d_C1_off
    d_TTR       = d_TTR     + f_d_C1_off
    d_TNFR      = d_TNFR    + f_d_C1_off
    
    d_C1        = d_C1      - f_d_C1
    d_TTR       = d_TTR     + f_d_C1
    d_TNFR      = d_TNFR    + f_d_C1
    
    # C1_off  <--> C1
    d_C1_off    = d_C1_off  - f_a_C1
    d_C1        = d_C1      + f_a_C1
    
    d_C1_off    = d_C1_off  + f_C1_A20
    d_C1        = d_C1      - f_C1_A20
    
    d_C1_off    = d_C1_off  + f_C1_off
    d_C1        = d_C1      - f_C1_off
    
    # C1 & C1_off internalization [treated as deg]
    d_C1        = d_C1      - f_i_C1
    d_C1_off    = d_C1_off  - f_i_C1_off
    
    # TNF-Dependent Activation
    
    # TNF degradation
    d_TNF       = d_TNF     - f_pd_m_tnf
    
    # TNF:TNFR binding
    d_TNFRtnf      = d_TNFRtnf    +   f_a_tnfrmtnf
    d_tnfrm     = d_tnfrm   - 3 * f_a_tnfrmtnf
    d_TNF       = d_TNF     -   f_a_tnfrmtnf
    
    d_TNFR      = d_TNFR    - f_a_TNFRtnf
    d_TNF       = d_TNF     - f_a_TNFRtnf
    d_TNFRtnf   = d_TNFRtnf + f_a_TNFRtnf
    
    d_TNFR      = d_TNFR    + f_d_TNFRtnf
    d_TNF       = d_TNF     + f_d_TNFRtnf
    d_TNFRtnf   = d_TNFRtnf - f_d_TNFRtnf
    
    # Complex I Generation  TTR + TNFRtnf <--> C1tnf_off; C1tnf
    d_C1tnf_off = d_C1tnf_off   + f_a_C1tnf_off
    d_TTR       = d_TTR         - f_a_C1tnf_off
    d_TNFRtnf   = d_TNFRtnf     - f_a_C1tnf_off
    
    d_C1tnf_off = d_C1tnf_off   - f_d_C1tnf_off
    d_TTR       = d_TTR         + f_d_C1tnf_off
    d_TNFRtnf   = d_TNFRtnf     + f_d_C1tnf_off
    
    d_C1tnf     = d_C1tnf       - f_d_C1tnf
    d_TTR       = d_TTR         + f_d_C1tnf
    d_TNFRtnf   = d_TNFRtnf     + f_d_C1tnf
    
    # C1tnf_off <--> C1tnf
    d_C1tnf_off = d_C1tnf_off   - f_a_C1tnf
    d_C1tnf     = d_C1tnf       + f_a_C1tnf
    
    d_C1tnf_off = d_C1tnf_off   + f_C1tnf_off
    d_C1tnf     = d_C1tnf       - f_C1tnf_off
    
    d_C1tnf_off = d_C1tnf_off   + f_C1tnf_A20
    d_C1tnf     = d_C1tnf       - f_C1tnf_A20
    
    # C1tnf & C1tnf_off internalization [treated as deg]
    d_C1tnf_off = d_C1tnf_off   - f_i_C1tnf_off
    d_C1tnf     = d_C1tnf       - f_i_C1tnf
    
    # C1tnf_off;C1tnf <--> C1_off;C1 + TNF
    d_C1tnf_off = d_C1tnf_off   - f_d_tnf_C1_off
    d_C1_off    = d_C1_off      + f_d_tnf_C1_off
    d_TNF       = d_TNF         + f_d_tnf_C1_off
    
    d_C1tnf_off = d_C1tnf_off   + f_a_tnf_C1_off
    d_C1_off    = d_C1_off      - f_a_tnf_C1_off
    d_TNF       = d_TNF         - f_a_tnf_C1_off
    
    d_C1tnf     = d_C1tnf       - f_d_tnf_C1
    d_C1        = d_C1          + f_d_tnf_C1
    d_TNF       = d_TNF         + f_d_tnf_C1
    
    d_C1tnf     = d_C1tnf       + f_a_tnf_C1
    d_C1        = d_C1          - f_a_tnf_C1
    d_TNF       = d_TNF         - f_a_tnf_C1
    
    # IKKK Regulation
    
    d_IKKK      = d_IKKK     + f_IKKK_on_C1
    d_IKKK_off  = d_IKKK_off - f_IKKK_on_C1
    
    d_IKKK      = d_IKKK     + f_IKKK_on_C1tnf
    d_IKKK_off  = d_IKKK_off - f_IKKK_on_C1tnf
    
    ## ADDED 
    # A20 mediated TRAF6 inactivation
    d_TRAF6    = d_TRAF6    + v.TP[15] * a20 * TRAF6s
    d_TRAF6s    = d_TRAF6s  - v.TP[15] * a20 * TRAF6s
    
    ## Save concentrations for next time step
    delta[1]      = d_IkBa
    delta[2]      = d_IkBan
    delta[3]      = d_IkBaNFkB
    delta[4]      = d_IkBaNFkBn
    delta[5]      = d_IkBat
    delta[6]      = d_IkBb
    delta[7]      = d_IkBbn
    delta[8]      = d_IkBbNFkB
    delta[9]      = d_IkBbNFkBn
    delta[10]     = d_IkBbt
    delta[11]     = d_IkBe
    delta[12]     = d_IkBen
    delta[13]     = d_IkBeNFkB
    delta[14]     = d_IkBeNFkBn
    delta[15]     = d_IkBet
    
    delta[16]     = d_NFkB
    delta[17]     = d_NFkBn
    # 
    delta[18]     = d_LPSo
    delta[19]     = d_LPS
    delta[20]     = d_LPSen
    delta[21]     = d_TLR4
    delta[22]     = d_TLR4en
    delta[23]     = d_TLR4LPS
    delta[24]     = d_TLR4LPSen
    delta[25]     = d_MyD88
    delta[26]     = d_MyD88s
    delta[27]     = d_TRIF
    delta[28]     = d_TRIFs
    delta[29]     = d_TRAF6
    delta[30]     = d_TRAF6s
    delta[31]     = d_IKKK_off
    delta[32]     = d_IKKK
    delta[33]     = d_IKK_off
    delta[34]     = d_IKK
    delta[35]     = d_IKK_i
    
    # TNF component 
    delta[36]     = d_TNFnas  
    delta[37]     = d_TNFmRNA 
    delta[38]     = d_TNFpro  
    delta[39]     = d_TNF  
    
    # TNFR part
    delta[40]     = d_tnfrm
    delta[41]     = d_TNFR
    delta[42]     = d_TNFRtnf
    delta[43]     = d_C1
    delta[44]     = d_C1_off  
    delta[45]     = d_C1tnf  
    delta[46]     = d_C1tnf_off   
    delta[47]     = d_TTR 
    
    # A20
    delta[48]     = d_a20 
    delta[49]     = d_a20t 
    
    delta
end