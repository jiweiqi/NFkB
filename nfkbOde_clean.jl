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
    flux_rsu_a         = v.NP[1] 
    flux_rsu_b         = v.NP[2] 
    flux_rsu_e         = v.NP[3] 
    flux_rsu_f         = v.TP[1] ; # tnf basal transcription [NEW]
    
    flux_rsr_an        = v.NP[4]   * (NFkBn^v.NP[7]) 
    flux_rsr_bn        = v.NP[5]   * (NFkBn^v.NP[8]) 
    flux_rsr_en        = v.NP[6]   * (NFkBn^v.NP[9]) 
    
    flux_rsr_fn        = v.TP[2]  * (NFkBn^v.TP[3]) / (NFkBn^v.TP[3] + v.TP[4]^v.TP[3]); # tnf transcription 
    flux_rp_f         = kpr * TNFnas; # tnf RNA processing
    
    
    flux_rd_a          = v.NP[13]  * IkBat	
    flux_rd_b          = v.NP[14]  * IkBbt	
    flux_rd_e          = v.NP[15]  * IkBet	
    flux_rd_f          = kdegm * TNFmRNA; # tnf mRNA degradation
    
    flux_ps_c_a        = ktl_a	* IkBat	
    flux_ps_c_b        = ktl_b	* IkBbt	
    flux_ps_c_e        = ktl_e	* IkBet	
    flux_ps_c_f        = ktl	* TNFmRNA;  # tnf protein synthsis 
    flux_pd_f          = kdeg_tnf    * TNFpro;  # tnf protein deg. (inside the cell)
    flux_sc_c_f        = ksec   * TNFpro;  # tnf protein secrection 
    
    
    flux_in_a          = v.NP[19]	* IkBa	
    flux_in_b          = v.NP[20]	* IkBb	
    flux_in_e          = v.NP[21]	* IkBe	
    flux_in_n          = v.NP[22]	* NFkB	
    flux_ex_a          = v.NP[23]	* IkBan	
    flux_ex_b          = v.NP[24]	* IkBbn	
    flux_ex_e          = v.NP[25]	* IkBen	
    flux_ex_n          = v.NP[26]	* NFkBn	
    flux_in_2an        = v.NP[27]	* IkBaNFkB	
    flux_in_2bn        = v.NP[28]	* IkBbNFkB	
    flux_in_2en        = v.NP[29]	* IkBeNFkB	
    flux_ex_2an        = v.NP[30]	* IkBaNFkBn	
    flux_ex_2bn        = v.NP[31]	* IkBbNFkBn	
    flux_ex_2en        = v.NP[32]	* IkBeNFkBn	
    
    flux_pd_c_a        = v.NP[33]	* IkBa	
    flux_pd_c_b        = v.NP[34]	* IkBb	
    flux_pd_c_e        = v.NP[35]	* IkBe	
    flux_pd_n_a        = v.NP[36]	* IkBan	
    flux_pd_n_b        = v.NP[37]	* IkBbn	
    flux_pd_n_e        = v.NP[38]	* IkBen	
    flux_pd_c_2an      = v.NP[39]	* IkBaNFkB	
    flux_pd_c_2bn      = v.NP[40]	* IkBbNFkB	
    flux_pd_c_2en      = v.NP[41]	* IkBeNFkB	
    flux_pd_n_2an      = v.NP[42]	* IkBaNFkBn	
    flux_pd_n_2bn      = v.NP[43]	* IkBbNFkBn	
    flux_pd_n_2en      = v.NP[44]	* IkBeNFkBn	
    
    flux_a_c_an        = v.NP[45]  * IkBa	* NFkB	
    flux_a_c_bn        = v.NP[46]  * IkBb	* NFkB	
    flux_a_c_en        = v.NP[47]  * IkBe	* NFkB	
    flux_a_n_an        = v.NP[48]	* IkBan	* NFkBn	
    flux_a_n_bn        = v.NP[49]	* IkBbn	* NFkBn	
    flux_a_n_en        = v.NP[50]  * IkBen	* NFkBn	
    flux_d_c_an        = v.NP[51]	* IkBaNFkB	
    flux_d_c_bn        = v.NP[52]	* IkBbNFkB	
    flux_d_c_en        = v.NP[53]	* IkBeNFkB	
    flux_d_n_an        = v.NP[54]	* IkBaNFkBn	
    flux_d_n_bn        = v.NP[55]	* IkBbNFkBn	
    flux_d_n_en        = v.NP[56]	* IkBeNFkBn	
    
    # IKK Mediated IkB Degradation [free & bound]
    flux_ph_c_a        = v.NP[57]  * IkBa * IKK_flux 
    flux_ph_c_b        = v.NP[58]	* IkBb * IKK_flux 
    flux_ph_c_e        = v.NP[59]	* IkBe * IKK_flux 
    
    flux_ph_c_an       = v.NP[60]  * IkBaNFkB * IKK_flux 
    flux_ph_c_bn       = v.NP[61]	* IkBbNFkB * IKK_flux 
    flux_ph_c_en       = v.NP[62]	* IkBeNFkB * IKK_flux 
    
    # A20 Fluxes 
    flux_rsu_a20       = v.NP[63] 
    flux_rsr_a20       = v.NP[64]  * (NFkBn^v.NP[65]) 
    flux_rd_a20        = v.NP[66]  * a20t 
    flux_ps_c_a20      = v.NP[67]  * NFkBn 
    flux_pd_c_a20      = v.NP[68]  * a20  
    
    if ( t .> v.NP[69])  # disables inducible A20 txn to match exp data
        flux_rsr_a20 *= 0.0
    end
    
    
    # -- IKK Activation Module & TLR4 module
    # LPS binding to the surface()
    flux_b_LPS         = v.IP[1]   * LPSo  
    
    # ligand binding
    flux_b_TLR4LPS     = v.IP[5]    * TLR4 * LPS 
    flux_ub_TLR4LPS    = v.IP[6]    * TLR4LPS 
    flux_b_TLR4LPSen   = v.IP[7]    * TLR4en * LPSen 
    flux_ub_TLR4LPSen  = v.IP[8]    * TLR4LPSen 
    
    # Activation modified on 0426
    flux_a_MyD88       = v.IP[18]    * (TLR4LPS)^v.IP[19] / ((TLR4LPS)^v.IP[19] + (v.IP[20])^v.IP[19])  * MyD88
    flux_i_MyD88       = v.IP[21]    * MyD88s 
    flux_a_TRIF        = v.IP[22]    * TLR4LPSen * TRIF
    flux_i_TRIF        = v.IP[23]   * TRIFs 
    
    # shuttling; modified on 0426
    flux_in_LPS        = kshut   * LPS 
    flux_out_LPS       = v.IP[3]   * LPSen
    flux_in_TLR4       = v.IP[12]   * TLR4
    flux_out_TLR4      = v.IP[13]   * TLR4en
    flux_in_TLR4LPS    = kshut_com   * TLR4LPS
    flux_out_TLR4LPS   = v.IP[15]   * TLR4LPSen
    
    # generation & degradation
    flux_g_TLR4        = ktlr4
    flux_d_TLR4        = v.IP[10]   * TLR4
    flux_d_TLR4en      = v.IP[11]   * TLR4en
    flux_d_LPSen       = v.IP[4]   * LPSen
    flux_d_TLR4LPSen   = v.IP[17]   * TLR4LPSen
    # modified 0515
    flux_d_TLR4LPS     = v.IP[16]   * TLR4LPS
    # -----------------------IKK module
    # TRAF6
    flux_a_TRAF6_MyD88s = v.IP[24]   * TRAF6 * MyD88s
    flux_a_TRAF6_TRIFs = v.IP[25]   * TRAF6 * TRIFs
    flux_i_TRAF6       = v.IP[26]   * TRAF6s
    
    # IKKK
    
    flux_IKKK_on       = v.IP[28]   * IKKK_off
    flux_IKKK_on_TRAF6s = v.IP[27]   * IKKK_off * TRAF6s
    flux_IKKK_off      = v.IP[29]   * IKKK
    
    # IKK
    flux_IKK_on        = v.IP[31]   * IKK_off
    flux_IKK_on_IKKK   = v.IP[30]   * IKK_off * IKKK
    flux_IKK_off       = v.IP[32]   * IKK
    flux_IKK_off_i     = v.IP[33]   * IKK
    flux_IKKi          = v.IP[34]   * IKK_i
    
    
    ## TNFR part 
    flux_pd_m_tnf      = v.AP[1]    * TNF ; # pd_m_tnf 45' half life of exogenous TNF
    flux_syn_tnfrm     = ktnfr;  # tnfrm synthesis [txn, tsl, localization]
    flux_pd_tnfrm      = v.AP[3]    * tnfrm;# tnfrm --> deg  -- 120' halflife
    flux_a_tnfrm       = v.AP[4]    * tnfrm;# 3tnfrm --> TNFR
    flux_d_TNFR        = v.AP[5]    * TNFR;# TNFR   --> 3tnfrm
    flux_i_TNFR        = v.AP[6]    * TNFR;# TNFR internalization -- 30' halflife
    flux_a_C1_off      = v.AP[7]    * TNFR * TTR;# TNFR + TTR --> C1_off
    flux_d_C1_off      = v.AP[8]    * C1_off;# C1_off --> TNFR + TTR
    flux_a_C1          = v.AP[9]    * C1_off;# C1_off --> C1
    flux_C1_off        = v.AP[10]   * C1;# C1     --> C1_off
    flux_C1_A20        = v.AP[11]   * C1 * a20;# C1     --> C1_off [A20 Mediated]
    flux_d_C1          = v.AP[12]   * C1; # UPDATE # C1     --> TNFR + TTR
    flux_i_C1_off      = v.AP[13]   * C1_off; # C1_off internalization
    flux_i_C1          = v.AP[14]   * C1; # C1 internalization
    flux_a_tnfrmtnf    = v.AP[15]   * tnfrm * TNF; # 3tnfrm + tnf --> TNFRtnf
    flux_a_TNFRtnf     = v.AP[16]   * TNFR * TNF; # TNFR + tnf --> TNFRtnf
    flux_d_TNFRtnf     = v.AP[17]   * TNFRtnf; # TNFRtnf   --> TNFR + tnf
    flux_i_TNFRtnf     = v.AP[18]   * TNFRtnf;  # TNFRtnf internalization -- 30' halflife
    
    flux_a_C1tnf_off   = v.AP[19]   * TNFRtnf * TTR; # TNFRtnf + TTR --> C1tnf_off
    flux_d_C1tnf_off   = v.AP[20]   * C1tnf_off; # C1tnf_off --> TNFRtnf + TTR
    flux_a_C1tnf       = v.AP[21]   * C1tnf_off; # C1tnf_off --> C1tnf
    flux_C1tnf_off     = v.AP[22]   * C1tnf; # C1tnf     --> C1tnf_off
    flux_C1tnf_A20     = v.AP[23]   * C1tnf * a20;  # C1tnf     --> C1tnf_off [A20 Mediated]
    flux_d_C1tnf       = v.AP[24]   * C1tnf; # updated  # C1tnf     --> TNFRtnf + TTR
    flux_i_C1tnf       = v.AP[25]   * C1tnf; # C1tnf internalization
    flux_i_C1tnf_off   = v.AP[26]   * C1tnf_off;  # C1tnf_off internalization
    flux_d_tnf_C1_off  = v.AP[27]   * C1tnf_off; # C1tnf_off --> C1_off + tnf
    flux_a_tnf_C1_off  = v.AP[28]   * C1_off * TNF; # C1_off + tnf --> C1tnf_off
    flux_d_tnf_C1      = v.AP[29]   * C1tnf; # C1tnf    --> C1 + tnf
    flux_a_tnf_C1      = v.AP[30]   * C1 * TNF; # C1 + tnf --> C1tnf
    flux_IKKK_on_C1    = v.AP[31]   * IKKK_off * C1; # IKKK_off --> IKKK [C1 mediated]?500
    flux_IKKK_on_C1tnf = v.AP[32]   * IKKK_off * C1tnf;# IKKK_off --> IKKK [C1tnf mediated]
    
    
    ## Set iteration changes [delta] to zero
    delta_IkBa        = 0.0
    delta_IkBan       = 0.0
    delta_IkBaNFkB    = 0.0
    delta_IkBaNFkBn   = 0.0
    delta_IkBat       = 0.0
    delta_IkBb        = 0.0
    delta_IkBbn       = 0.0
    delta_IkBbNFkB    = 0.0
    delta_IkBbNFkBn   = 0.0
    delta_IkBbt       = 0.0
    delta_IkBe        = 0.0
    delta_IkBen       = 0.0
    delta_IkBeNFkB    = 0.0
    delta_IkBeNFkBn   = 0.0
    delta_IkBet       = 0.0
    delta_NFkB        = 0.0
    delta_NFkBn       = 0.0
    # TLR4 module
    delta_LPS         = 0.0
    delta_TLR4        = 0.0
    delta_TLR4LPS     = 0.0
    delta_MyD88       = 0.0
    delta_MyD88s      = 0.0
    delta_TRIF        = 0.0
    delta_TRIFs       = 0.0
    delta_LPSen       = 0.0
    delta_TLR4LPSen   = 0.0
    # IKK module
    delta_IKKK        = 0.0
    delta_IKKK_off    = 0.0
    delta_IKK         = 0.0
    delta_IKK_off     = 0.0
    delta_IKK_i       = 0.0
    # TLR4 continue()
    delta_TRAF6       = 0.0
    delta_TRAF6s      = 0.0
    delta_TLR4en      = 0.0
    
    
    delta_LPSo           = 0.0
    
    # TNF component 
    delta_TNFnas  = 0.0 
    delta_TNFmRNA = 0.0 
    delta_TNFpro  = 0.0 
    delta_TNF  = 0.0 
    
    # TNFR part 
    delta_tnfrm       = 0.0
    delta_TNFR        = 0.0
    delta_TNFRtnf     = 0.0
    delta_C1          = 0.0
    delta_C1_off      = 0.0
    delta_C1tnf       = 0.0
    delta_C1tnf_off   = 0.0
    delta_TTR         = 0.0
    
    # A20 NEW
    delta_a20         = 0.0
    delta_a20t        = 0.0
    
    
    ## Add Fluxes to component concentrations; then Save
    
    # IkB Transcrv.IPtion
    delta_IkBat     = delta_IkBat + flux_rsu_a
    delta_IkBbt     = delta_IkBbt + flux_rsu_b
    delta_IkBet     = delta_IkBet + flux_rsu_e
    delta_IkBat     = delta_IkBat + flux_rsr_an
    delta_IkBbt     = delta_IkBbt + flux_rsr_bn
    delta_IkBet     = delta_IkBet + flux_rsr_en
    delta_TNFnas    = delta_TNFnas + flux_rsu_f; # NEW
    delta_TNFnas    = delta_TNFnas + flux_rsr_fn; # NEW
    delta_TNFnas    = delta_TNFnas - flux_rp_f; # NEW
    delta_TNFmRNA   = delta_TNFmRNA + flux_rp_f; # NEW
    
    
    # IkB Transcrv.IPt Degradation
    delta_IkBbt     = delta_IkBbt - flux_rd_b
    delta_IkBet     = delta_IkBet - flux_rd_e
    delta_IkBat     = delta_IkBat - flux_rd_a
    delta_TNFmRNA   = delta_TNFmRNA  - flux_rd_f;   # tnft degradation
    
    # IkB Translation
    delta_IkBat     = delta_IkBat - flux_ps_c_a
    delta_IkBat     = delta_IkBat + flux_ps_c_a
    delta_IkBa      = delta_IkBa  + flux_ps_c_a
    
    delta_IkBbt     = delta_IkBbt - flux_ps_c_b
    delta_IkBbt     = delta_IkBbt + flux_ps_c_b
    delta_IkBb      = delta_IkBb  + flux_ps_c_b
    
    delta_IkBet     = delta_IkBet - flux_ps_c_e
    delta_IkBet     = delta_IkBet + flux_ps_c_e
    delta_IkBe      = delta_IkBe  + flux_ps_c_e
    
    delta_TNFpro    = delta_TNFpro + flux_ps_c_f; 
    
    # TNF secretion + protein degradation
    delta_TNFpro     = delta_TNFpro - flux_sc_c_f; 
    delta_TNF        = delta_TNF    + flux_sc_c_f; 
    delta_TNFpro     = delta_TNFpro - flux_pd_f; 
    
    # IkB:NFkB Shuttling [Free & Bound]
    delta_IkBa      = delta_IkBa  - flux_in_a
    delta_IkBan     = delta_IkBan + flux_in_a
    
    delta_IkBb      = delta_IkBb  - flux_in_b
    delta_IkBbn     = delta_IkBbn + flux_in_b
    
    delta_IkBe      = delta_IkBe  - flux_in_e
    delta_IkBen     = delta_IkBen + flux_in_e
    
    delta_NFkB      = delta_NFkB  - flux_in_n
    delta_NFkBn     = delta_NFkBn + flux_in_n
    
    delta_IkBan     = delta_IkBan - flux_ex_a
    delta_IkBa      = delta_IkBa  + flux_ex_a
    
    delta_IkBbn     = delta_IkBbn - flux_ex_b
    delta_IkBb      = delta_IkBb  + flux_ex_b
    
    delta_IkBen     = delta_IkBen - flux_ex_e
    delta_IkBe      = delta_IkBe  + flux_ex_e
    
    delta_NFkBn     = delta_NFkBn - flux_ex_n
    delta_NFkB      = delta_NFkB  + flux_ex_n
    
    delta_IkBaNFkB  = delta_IkBaNFkB  - flux_in_2an
    delta_IkBaNFkBn = delta_IkBaNFkBn + flux_in_2an
    
    delta_IkBbNFkB  = delta_IkBbNFkB  - flux_in_2bn
    delta_IkBbNFkBn = delta_IkBbNFkBn + flux_in_2bn
    
    delta_IkBeNFkB  = delta_IkBeNFkB  - flux_in_2en
    delta_IkBeNFkBn = delta_IkBeNFkBn + flux_in_2en
    
    delta_IkBaNFkBn = delta_IkBaNFkBn - flux_ex_2an
    delta_IkBaNFkB  = delta_IkBaNFkB  + flux_ex_2an
    
    delta_IkBbNFkBn = delta_IkBbNFkBn - flux_ex_2bn
    delta_IkBbNFkB  = delta_IkBbNFkB  + flux_ex_2bn
    
    delta_IkBeNFkBn = delta_IkBeNFkBn - flux_ex_2en
    delta_IkBeNFkB  = delta_IkBeNFkB  + flux_ex_2en
    
    # IkB:NFkB Association [Cytoplasm & Nucleus]
    delta_IkBa      = delta_IkBa - flux_a_c_an
    delta_NFkB      = delta_NFkB - flux_a_c_an
    delta_IkBaNFkB  = delta_IkBaNFkB + flux_a_c_an
    
    delta_IkBb      = delta_IkBb - flux_a_c_bn
    delta_NFkB      = delta_NFkB - flux_a_c_bn
    delta_IkBbNFkB  = delta_IkBbNFkB + flux_a_c_bn
    
    delta_IkBe      = delta_IkBe - flux_a_c_en
    delta_NFkB      = delta_NFkB - flux_a_c_en
    delta_IkBeNFkB  = delta_IkBeNFkB + flux_a_c_en
    
    delta_IkBan     = delta_IkBan - flux_a_n_an
    delta_NFkBn     = delta_NFkBn - flux_a_n_an
    delta_IkBaNFkBn = delta_IkBaNFkBn + flux_a_n_an
    
    delta_IkBbn     = delta_IkBbn - flux_a_n_bn
    delta_NFkBn     = delta_NFkBn - flux_a_n_bn
    delta_IkBbNFkBn = delta_IkBbNFkBn + flux_a_n_bn
    
    delta_IkBen     = delta_IkBen - flux_a_n_en
    delta_NFkBn     = delta_NFkBn - flux_a_n_en
    delta_IkBeNFkBn = delta_IkBeNFkBn + flux_a_n_en
    
    # IkB:NFkB Dissociation [Cytoplasm & Nucleus]
    delta_IkBaNFkB  = delta_IkBaNFkB - flux_d_c_an
    delta_IkBa      = delta_IkBa + flux_d_c_an
    delta_NFkB      = delta_NFkB + flux_d_c_an
    
    delta_IkBbNFkB  = delta_IkBbNFkB - flux_d_c_bn
    delta_IkBb      = delta_IkBb + flux_d_c_bn
    delta_NFkB      = delta_NFkB + flux_d_c_bn
    
    delta_IkBeNFkB  = delta_IkBeNFkB - flux_d_c_en
    delta_IkBe      = delta_IkBe + flux_d_c_en
    delta_NFkB      = delta_NFkB + flux_d_c_en
    
    delta_IkBaNFkBn = delta_IkBaNFkBn - flux_d_n_an
    delta_IkBan     = delta_IkBan + flux_d_n_an
    delta_NFkBn     = delta_NFkBn + flux_d_n_an
    
    delta_IkBbNFkBn = delta_IkBbNFkBn - flux_d_n_bn
    delta_IkBbn     = delta_IkBbn + flux_d_n_bn
    delta_NFkBn     = delta_NFkBn + flux_d_n_bn
    
    delta_IkBeNFkBn = delta_IkBeNFkBn - flux_d_n_en
    delta_IkBen     = delta_IkBen + flux_d_n_en
    delta_NFkBn     = delta_NFkBn + flux_d_n_en
    
    # Free IkB Degradation [Cytoplasm & Nucleus]
    delta_IkBa      = delta_IkBa  - flux_pd_c_a
    delta_IkBb      = delta_IkBb  - flux_pd_c_b
    delta_IkBe      = delta_IkBe  - flux_pd_c_e
    delta_IkBan     = delta_IkBan - flux_pd_n_a
    delta_IkBbn     = delta_IkBbn - flux_pd_n_b
    delta_IkBen     = delta_IkBen - flux_pd_n_e
    
    # IkB:NFkB Degradation [Cytoplasm & Nucleus]
    delta_IkBaNFkB  = delta_IkBaNFkB - flux_pd_c_2an
    delta_NFkB      = delta_NFkB + flux_pd_c_2an
    
    delta_IkBbNFkB  = delta_IkBbNFkB - flux_pd_c_2bn
    delta_NFkB      = delta_NFkB + flux_pd_c_2bn
    
    delta_IkBeNFkB  = delta_IkBeNFkB - flux_pd_c_2en
    delta_NFkB      = delta_NFkB + flux_pd_c_2en
    
    delta_IkBaNFkBn = delta_IkBaNFkBn - flux_pd_n_2an
    delta_NFkBn     = delta_NFkBn + flux_pd_n_2an
    
    delta_IkBbNFkBn = delta_IkBbNFkBn - flux_pd_n_2bn
    delta_NFkBn     = delta_NFkBn + flux_pd_n_2bn
    
    delta_IkBeNFkBn = delta_IkBeNFkBn - flux_pd_n_2en
    delta_NFkBn     = delta_NFkBn + flux_pd_n_2en
    
    # IKK Mediated IkB Degradation
    delta_IkBa      = delta_IkBa - flux_ph_c_a
    delta_IkBb      = delta_IkBb - flux_ph_c_b
    delta_IkBe      = delta_IkBe - flux_ph_c_e
    
    # IKK Mediated IkB:NFkB Degradation
    delta_IkBaNFkB  = delta_IkBaNFkB - flux_ph_c_an
    delta_NFkB      = delta_NFkB + flux_ph_c_an
    
    delta_IkBbNFkB  = delta_IkBbNFkB - flux_ph_c_bn
    delta_NFkB      = delta_NFkB + flux_ph_c_bn
    
    delta_IkBeNFkB  = delta_IkBeNFkB - flux_ph_c_en
    delta_NFkB      = delta_NFkB + flux_ph_c_en
    
    
    # Upstream Pathway  [LPS] --> IKK activation
    
    # Ligand binding
    delta_LPSo      = delta_LPSo - flux_b_LPS
    delta_LPS       = delta_LPS  + flux_b_LPS
    
    delta_LPS       = delta_LPS       -   flux_b_TLR4LPS
    delta_TLR4       = delta_TLR4       -   flux_b_TLR4LPS
    delta_TLR4LPS   = delta_TLR4LPS   +   flux_b_TLR4LPS
    delta_LPS       = delta_LPS       +   flux_ub_TLR4LPS
    delta_TLR4       = delta_TLR4     +   flux_ub_TLR4LPS
    delta_TLR4LPS   = delta_TLR4LPS   -   flux_ub_TLR4LPS
    
    delta_LPSen     = delta_LPSen     -   flux_b_TLR4LPSen
    delta_TLR4en     = delta_TLR4en     -   flux_b_TLR4LPSen
    delta_TLR4LPSen = delta_TLR4LPSen +   flux_b_TLR4LPSen
    delta_LPSen     = delta_LPSen     +   flux_ub_TLR4LPSen
    delta_TLR4en     = delta_TLR4en   +   flux_ub_TLR4LPSen
    delta_TLR4LPSen = delta_TLR4LPSen -   flux_ub_TLR4LPSen
    
    # Activation
    delta_MyD88     = delta_MyD88     -   flux_a_MyD88
    delta_MyD88s    = delta_MyD88s    +   flux_a_MyD88
    delta_MyD88     = delta_MyD88     +   flux_i_MyD88
    delta_MyD88s    = delta_MyD88s    -   flux_i_MyD88
    
    delta_TRIF      = delta_TRIF      -   flux_a_TRIF
    delta_TRIFs     = delta_TRIFs     +   flux_a_TRIF
    delta_TRIF      = delta_TRIF      +   flux_i_TRIF
    delta_TRIFs     = delta_TRIFs     -   flux_i_TRIF
    
    # shuttling
    delta_LPS       = delta_LPS       -   flux_in_LPS
    delta_LPSen     = delta_LPSen     +   flux_in_LPS
    # delta_LPSo       = delta_LPSo     +   flux_out_LPS
    delta_LPS       = delta_LPS       +   flux_out_LPS
    delta_LPSen     = delta_LPSen     -   flux_out_LPS
    
    delta_TLR4      = delta_TLR4      -   flux_in_TLR4
    delta_TLR4en    = delta_TLR4en    +   flux_in_TLR4
    delta_TLR4      = delta_TLR4      +   flux_out_TLR4
    delta_TLR4en    = delta_TLR4en    -   flux_out_TLR4
    
    delta_TLR4LPS      = delta_TLR4LPS      -   flux_in_TLR4LPS
    delta_TLR4LPSen    = delta_TLR4LPSen    +   flux_in_TLR4LPS
    delta_TLR4LPS      = delta_TLR4LPS      +   flux_out_TLR4LPS
    delta_TLR4LPSen     = delta_TLR4LPSen   -   flux_out_TLR4LPS
    
    
    # generation & degradation
    delta_TLR4      = delta_TLR4 + flux_g_TLR4
    delta_TLR4      = delta_TLR4 - flux_d_TLR4
    delta_TLR4en     = delta_TLR4en - flux_d_TLR4en
    delta_LPSen      = delta_LPSen - flux_d_LPSen
    delta_TLR4LPSen      = delta_TLR4LPSen - flux_d_TLR4LPSen
    # -- modified 0515
    delta_TLR4LPS      = delta_TLR4LPS - flux_d_TLR4LPS
    
    
    # TRAF6
    delta_TRAF6       = delta_TRAF6  - flux_a_TRAF6_MyD88s
    delta_TRAF6s      = delta_TRAF6s + flux_a_TRAF6_MyD88s
    delta_TRAF6       = delta_TRAF6  - flux_a_TRAF6_TRIFs
    delta_TRAF6s      = delta_TRAF6s + flux_a_TRAF6_TRIFs
    delta_TRAF6       = delta_TRAF6  + flux_i_TRAF6
    delta_TRAF6s      = delta_TRAF6s - flux_i_TRAF6
    
    # IKKK Regulation
    delta_IKKK      = delta_IKKK     + flux_IKKK_on
    delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on
    
    delta_IKKK      = delta_IKKK     + flux_IKKK_on_TRAF6s
    delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on_TRAF6s
    
    delta_IKKK_off  = delta_IKKK_off + flux_IKKK_off
    delta_IKKK      = delta_IKKK     - flux_IKKK_off
    
    # IKK Regulation
    delta_IKK       = delta_IKK     + flux_IKK_on
    delta_IKK_off   = delta_IKK_off - flux_IKK_on
    
    delta_IKK       = delta_IKK     + flux_IKK_on_IKKK
    delta_IKK_off   = delta_IKK_off - flux_IKK_on_IKKK
    
    delta_IKK       = delta_IKK     - flux_IKK_off
    delta_IKK_off   = delta_IKK_off + flux_IKK_off
    
    delta_IKK       = delta_IKK     - flux_IKK_off_i
    delta_IKK_i     = delta_IKK_i   + flux_IKK_off_i
    
    delta_IKK_off   = delta_IKK_off + flux_IKKi
    delta_IKK_i     = delta_IKK_i   - flux_IKKi
    
    
    
    
    
    
    ## TNR module [NEW] + A20 [NEW]
    # A20 Expression Transcrviption
    delta_a20t      = delta_a20t + flux_rsu_a20
    delta_a20t      = delta_a20t + flux_rsr_a20
    delta_a20t      = delta_a20t - flux_rd_a20
    
    # Translation
    delta_a20       = delta_a20  + flux_ps_c_a20
    
    # Degradation
    delta_a20       = delta_a20 - flux_pd_c_a20
    
    # Upstream Pathway  [TNF] --> IKK activation
    
    # TNF-independent Activation
    
    # tnfrm trimerization 3tnfrm<--> TNFR
    delta_tnfrm     = delta_tnfrm - 3 * flux_a_tnfrm
    delta_TNFR      = delta_TNFR +   flux_a_tnfrm
    
    delta_tnfrm     = delta_tnfrm + 3 * flux_d_TNFR
    delta_TNFR      = delta_TNFR -   flux_d_TNFR
    
    # TNF Receptor Metabolism
    delta_TNFR     = delta_TNFR     - flux_i_TNFR
    delta_TNFRtnf  = delta_TNFRtnf  - flux_i_TNFRtnf
    
    delta_tnfrm    = delta_tnfrm    - flux_pd_tnfrm
    delta_tnfrm    = delta_tnfrm    + flux_syn_tnfrm
    
    # Complex I Generation  TTR + TNFR <--> C1_off; C1
    delta_C1_off    = delta_C1_off  + flux_a_C1_off
    delta_TTR       = delta_TTR     - flux_a_C1_off
    delta_TNFR      = delta_TNFR    - flux_a_C1_off
    
    delta_C1_off    = delta_C1_off  - flux_d_C1_off
    delta_TTR       = delta_TTR     + flux_d_C1_off
    delta_TNFR      = delta_TNFR    + flux_d_C1_off
    
    delta_C1        = delta_C1      - flux_d_C1
    delta_TTR       = delta_TTR     + flux_d_C1
    delta_TNFR      = delta_TNFR    + flux_d_C1
    
    # C1_off  <--> C1
    delta_C1_off    = delta_C1_off  - flux_a_C1
    delta_C1        = delta_C1      + flux_a_C1
    
    delta_C1_off    = delta_C1_off  + flux_C1_A20
    delta_C1        = delta_C1      - flux_C1_A20
    
    delta_C1_off    = delta_C1_off  + flux_C1_off
    delta_C1        = delta_C1      - flux_C1_off
    
    # C1 & C1_off internalization [treated as deg]
    delta_C1        = delta_C1      - flux_i_C1
    delta_C1_off    = delta_C1_off  - flux_i_C1_off
    
    # TNF-Dependent Activation
    
    # TNF degradation
    delta_TNF       = delta_TNF     - flux_pd_m_tnf
    
    # TNF:TNFR binding
    delta_TNFRtnf      = delta_TNFRtnf    +   flux_a_tnfrmtnf
    delta_tnfrm     = delta_tnfrm   - 3 * flux_a_tnfrmtnf
    delta_TNF       = delta_TNF     -   flux_a_tnfrmtnf
    
    delta_TNFR      = delta_TNFR    - flux_a_TNFRtnf
    delta_TNF       = delta_TNF     - flux_a_TNFRtnf
    delta_TNFRtnf   = delta_TNFRtnf + flux_a_TNFRtnf
    
    delta_TNFR      = delta_TNFR    + flux_d_TNFRtnf
    delta_TNF       = delta_TNF     + flux_d_TNFRtnf
    delta_TNFRtnf   = delta_TNFRtnf - flux_d_TNFRtnf
    
    # Complex I Generation  TTR + TNFRtnf <--> C1tnf_off; C1tnf
    delta_C1tnf_off = delta_C1tnf_off   + flux_a_C1tnf_off
    delta_TTR       = delta_TTR         - flux_a_C1tnf_off
    delta_TNFRtnf   = delta_TNFRtnf     - flux_a_C1tnf_off
    
    delta_C1tnf_off = delta_C1tnf_off   - flux_d_C1tnf_off
    delta_TTR       = delta_TTR         + flux_d_C1tnf_off
    delta_TNFRtnf   = delta_TNFRtnf     + flux_d_C1tnf_off
    
    delta_C1tnf     = delta_C1tnf       - flux_d_C1tnf
    delta_TTR       = delta_TTR         + flux_d_C1tnf
    delta_TNFRtnf   = delta_TNFRtnf     + flux_d_C1tnf
    
    # C1tnf_off <--> C1tnf
    delta_C1tnf_off = delta_C1tnf_off   - flux_a_C1tnf
    delta_C1tnf     = delta_C1tnf       + flux_a_C1tnf
    
    delta_C1tnf_off = delta_C1tnf_off   + flux_C1tnf_off
    delta_C1tnf     = delta_C1tnf       - flux_C1tnf_off
    
    delta_C1tnf_off = delta_C1tnf_off   + flux_C1tnf_A20
    delta_C1tnf     = delta_C1tnf       - flux_C1tnf_A20
    
    # C1tnf & C1tnf_off internalization [treated as deg]
    delta_C1tnf_off = delta_C1tnf_off   - flux_i_C1tnf_off
    delta_C1tnf     = delta_C1tnf       - flux_i_C1tnf
    
    # C1tnf_off;C1tnf <--> C1_off;C1 + TNF
    delta_C1tnf_off = delta_C1tnf_off   - flux_d_tnf_C1_off
    delta_C1_off    = delta_C1_off      + flux_d_tnf_C1_off
    delta_TNF       = delta_TNF         + flux_d_tnf_C1_off
    
    delta_C1tnf_off = delta_C1tnf_off   + flux_a_tnf_C1_off
    delta_C1_off    = delta_C1_off      - flux_a_tnf_C1_off
    delta_TNF       = delta_TNF         - flux_a_tnf_C1_off
    
    delta_C1tnf     = delta_C1tnf       - flux_d_tnf_C1
    delta_C1        = delta_C1          + flux_d_tnf_C1
    delta_TNF       = delta_TNF         + flux_d_tnf_C1
    
    delta_C1tnf     = delta_C1tnf       + flux_a_tnf_C1
    delta_C1        = delta_C1          - flux_a_tnf_C1
    delta_TNF       = delta_TNF         - flux_a_tnf_C1
    
    # IKKK Regulation
    
    delta_IKKK      = delta_IKKK     + flux_IKKK_on_C1
    delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on_C1
    
    delta_IKKK      = delta_IKKK     + flux_IKKK_on_C1tnf
    delta_IKKK_off  = delta_IKKK_off - flux_IKKK_on_C1tnf
    
    
    
    ## ADDED 
    # A20 mediated TRAF6 inactivation
    delta_TRAF6    = delta_TRAF6    + v.TP[15] * a20 * TRAF6s
    delta_TRAF6s    = delta_TRAF6s  - v.TP[15] * a20 * TRAF6s
    
    
    ## Save concentrations for next time step
    delta[1]      = delta_IkBa
    delta[2]      = delta_IkBan
    delta[3]      = delta_IkBaNFkB
    delta[4]      = delta_IkBaNFkBn
    delta[5]      = delta_IkBat
    delta[6]      = delta_IkBb
    delta[7]      = delta_IkBbn
    delta[8]      = delta_IkBbNFkB
    delta[9]      = delta_IkBbNFkBn
    delta[10]     = delta_IkBbt
    delta[11]     = delta_IkBe
    delta[12]     = delta_IkBen
    delta[13]     = delta_IkBeNFkB
    delta[14]     = delta_IkBeNFkBn
    delta[15]     = delta_IkBet
    
    
    
    delta[16]     = delta_NFkB
    delta[17]     = delta_NFkBn
    # 
    delta[18]     = delta_LPSo
    delta[19]     = delta_LPS
    delta[20]     = delta_LPSen
    delta[21]     = delta_TLR4
    delta[22]     = delta_TLR4en
    delta[23]     = delta_TLR4LPS
    delta[24]     = delta_TLR4LPSen
    delta[25]     = delta_MyD88
    delta[26]     = delta_MyD88s
    delta[27]     = delta_TRIF
    delta[28]     = delta_TRIFs
    delta[29]     = delta_TRAF6
    delta[30]     = delta_TRAF6s
    delta[31]     = delta_IKKK_off
    delta[32]     = delta_IKKK
    delta[33]     = delta_IKK_off
    delta[34]     = delta_IKK
    delta[35]     = delta_IKK_i
    
    
    # TNF component 
    delta[36]     = delta_TNFnas  
    delta[37]     = delta_TNFmRNA 
    delta[38]     = delta_TNFpro  
    delta[39]     = delta_TNF  
    
    # TNFR part
    delta[40]     = delta_tnfrm
    delta[41]     = delta_TNFR
    delta[42]     = delta_TNFRtnf
    delta[43]     = delta_C1
    delta[44]     = delta_C1_off  
    delta[45]     = delta_C1tnf  
    delta[46]     = delta_C1tnf_off   
    delta[47]     = delta_TTR 
    
    # A20
    delta[48]     = delta_a20 
    delta[49]     = delta_a20t 
    
    delta
    # # if (v.flag_noTnfFeedback) && (t>0)
    # #     v.TP[8]=kdeg_tnf
    # # end 
end