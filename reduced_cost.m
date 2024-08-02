%% Hayden Gallo
% Bucci Lab
% Cobratoolbox Modeling
% 7/29/24
% Determining reduced cost of different medias to obtain growth from F.
% prau and P. copri


%% initialize cobratoolbox

initCobraToolbox


%% Load different models 

f_prau_model = readCbModel('/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/Faecalibacterium_prausnitzii_ERR1022327.mat');

p_copri_model = readCbModel('/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/Prevotella_copri_ERR1022397.mat');
p_copri_model_alt = readCbModel('/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/Prevotella_copri_CB7_DSM_18205.mat');

combined_model = readCbModel('/Users/haydengallo/cobratoolbox/Results_p_copri_f_prau/pairedModel_Prevotella_copri_ERR1022397_Faecalibacterium_prausnitzii_ERR1022484.mat');

%% Diets/Medias

pyg_media = {'EX_mobd(e)','-3.07E-04','1000';'EX_ac(e)','-0.0','1000';'EX_n2(e)','-227.8571429','1000';'EX_na1(e)','-20.51163798','1000';'EX_cl(e)','-5.941065976','1000';'EX_ca2(e)','-0.173941043','1000';'EX_fe2(e)','-0.016053362','1000';'EX_mg2(e)','-0.474477191','1000';'EX_k(e)','-35.39748582','1000';'EX_so4(e)','-0.710983412','1000';'EX_pi(e)','-18.29826648','1000';'EX_ala_L(e)','-16.89227108','1000';'EX_arg_L(e)','-5.338568575','1000';'EX_asn_L(e)','-1.286718791','1000';'EX_asp_L(e)','-10.81893313','1000';'EX_Lcystin(e)','-0.187272152','1000';'EX_glu_L(e)','-18.56754922','1000';'EX_gln_L(e)','-0.205274178','1000';'EX_gly(e)','-20.3151851','1000';'EX_his_L(e)','-3.319218598','1000';'EX_ile_L(e)','-8.195159139','1000';'EX_leu_L(e)','-11.2826377','1000';'EX_lys_l(e)','-9.884397018','1000';'EX_met_L(e)','-2.144657123','1000';'EX_phe_L(e)','-6.083829725','1000';'EX_pro_L(e)','-11.89938505','1000';'EX_ser_L(e)','-4.424652451','1000';'EX_thr_L(e)','-3.567830759','1000';'EX_trp_L(e)','-0.685504997','1000';'EX_tyr_L(e)','-1.683306566','1000';'EX_val_L(e)','-10.37149589','1000';'EX_glc_D(e)','-27.7537245498346','1000';'EX_hco3(e)','-1.190379826','1000';'EX_phyQ(e)','-0.002172408','1000';'EX_etoh(e)','-3.237913066','1000';'EX_pheme(e)','-0.0076693','1000';'EX_oh1(e)','-0.099987502','1000';'EX_cys_L(e)','-2.846894039','1000';'EX_M02144(e)','-2.846894039','1000';'EX_h2o(e)','-55013.2623','1000';'EX_h(e)','-6.30957E-05','1000'};
pyg_test = {'EX_n2(e)','-227.8571429','1000';'EX_na1(e)','-20.51163798','1000';'EX_cl(e)','-5.941065976','1000';'EX_ca2(e)','-0.173941043','1000';'EX_fe2(e)','-0.016053362','1000';'EX_mg2(e)','-0.474477191','1000';'EX_k(e)','-35.39748582','1000';'EX_so4(e)','-0.710983412','1000';'EX_pi(e)','-18.29826648','1000';'EX_ala_L(e)','-16.89227108','1000';'EX_arg_L(e)','-5.338568575','1000';'EX_asn_L(e)','-1.286718791','1000';'EX_asp_L(e)','-10.81893313','1000';'EX_Lcystin(e)','-0.187272152','1000';'EX_glu_L(e)','-18.56754922','1000';'EX_gln_L(e)','-0.205274178','1000';'EX_gly(e)','-20.3151851','1000';'EX_his_L(e)','-3.319218598','1000';'EX_ile_L(e)','-8.195159139','1000';'EX_leu_L(e)','-11.2826377','1000';'EX_lys_l(e)','-9.884397018','1000';'EX_met_L(e)','-2.144657123','1000';'EX_phe_L(e)','-6.083829725','1000';'EX_pro_L(e)','-11.89938505','1000';'EX_ser_L(e)','-4.424652451','1000';'EX_thr_L(e)','-3.567830759','1000';'EX_trp_L(e)','-0.685504997','1000';'EX_tyr_L(e)','-1.683306566','1000';'EX_val_L(e)','-10.37149589','1000';'EX_glc_D(e)','-27.7537245498346','1000';'EX_hco3(e)','-1.190379826','1000';'EX_phyQ(e)','-0.002172408','1000';'EX_etoh(e)','-3.237913066','1000';'EX_pheme(e)','-0.0076693','1000';'EX_oh1(e)','-0.099987502','1000';'EX_cys_L(e)','-2.846894039','1000';'EX_M02144(e)','-2.846894039','1000';'EX_h2o(e)','-55013.2623','1000';'EX_h(e)','-6.30957E-05','1000'};
ABB_media={'EX_pi(e)','-3.777186','1000';'EX_k(e)','-12.084479','1000';'EX_ac(e)','-0.0','1000';'EX_arg_L(e)','-9.357036999999998','1000';'EX_arg_D(e)','-9.357036999999998','1000';'EX_btn(e)','-3.86e-05','1000';'EX_cl(e)','-90.79522399999998','1000';'EX_cbl1(e)','-10e-4','1000';'EX_co2(e)','-100.0','1000';'EX_cys_L(e)','-3.8169380000000004','1000';'EX_h(e)','-3.172326','1000';'EX_dtt(e)','-6.482982','1000';'EX_so4(e)','-3.9465080000000006','1000';'EX_h2o(e)','-77713.01693','1000';'EX_na1(e)','-131.369616','1000';'EX_glc_D(e)','-5.550622','1000';'EX_h2(e)','-50.0','1000';'EX_mg2(e)','-3.393075','1000';'EX_mqn8(e)','-0.0011200000000000001','1000';'EX_mn2(e)','-0.01','1000';'EX_n2(e)','-850.0','1000';'EX_hco3(e)','-4.761338','1000';'EX_pyr(e)','-9.087605','1000';'EX_succ(e)','-1.850892','1000';'EX_nh4(e)','-3.2222280000000003','1000';'EX_ppi(e)','-0.6709520000000001','1000';'EX_fe3(e)','-0.731952','1000';'EX_strch1(e)','-0.504602','1000';'EX_thiog(e)','-4.382121','1000';'EX_thm(e)','-0.0029000000000000002','1000';'EX_ala_D(e)','-8.889688999999999','1000';'EX_ala_L(e)','-8.889688999999999','1000';'EX_asp_L(e)','-9.718608','1000';'EX_gln_L(e)','-0.314768','1000';'EX_gln_D(e)','-0.314768','1000';'EX_glu_D(e)','-12.081918','1000';'EX_glu_L(e)','-12.081918','1000';'EX_gly(e)','-7.273324000000001','1000';'EX_his_L(e)','-1.7208890000000001','1000';'EX_ile_L(e)','-5.01633','1000';'EX_leu_L(e)','-7.433019','1000';'EX_lys_L(e)','-5.376564999999999','1000';'EX_met_L(e)','-0.911472','1000';'EX_phe_L(e)','-4.104351','1000';'EX_pro_L(e)','-3.9954629999999995','1000';'EX_ser_L(e)','-3.3494800000000002','1000';'EX_asn_L(e)','-2.0073269999999996','1000';'EX_trp_L(e)','-0.328062','1000';'EX_tyr_L(e)','-1.6115650000000001','1000';'EX_val_L(e)','-5.778952','1000';'EX_ca2(e)','-0.242261','1000';'EX_fe2(e)','-0.061','1000';'EX_thr_L(e)','-0.940228','1000';'EX_cd2(e)','-9.340000000000001e-05','1000';'EX_cobalt2(e)','-0.01','1000';'EX_cu2(e)','-0.01','1000';'EX_ni2(e)','-0.000763','1000';'EX_zn2(e)','-0.395929','1000';'EX_ade(e)','-0.457929','1000';'EX_gua(e)','-0.428442','1000';'EX_csn(e)','-0.274698','1000';'EX_ura(e)','-0.394065','1000';'EX_man(e)','-1.779569','1000';'EX_pnto_R(e)','-0.0033799999999999998','1000';'EX_ascb_L(e)','-0.00011899999999999999','1000';'EX_nac(e)','-0.0229','1000';'EX_pydx(e)','-0.0006209999999999999','1000';'EX_chol(e)','-0.021500000000000002','1000';'EX_adocbl(e)','-3.62e-09','1000';'EX_mobd(e)','-0.00043','1000';'EX_26dap_M(e)','-10e-4','1000';'EX_h2s(e)','-10e-4','1000';'EX_ribflv(e)','-10e-4','1000';'EX_thymd(e)','-10e-4','1000';'EX_fol(e)','-10e-4','1000';'EX_12dgr180(e)','-10e-4','1000';'EX_hxan(e)','-10e-4','1000';'EX_ddca(e)','-10e-4','1000';'EX_ttdca(e)','-10e-4','1000';'EX_pheme(e)','-10e-4','1000';'EX_q8(e)','-10e-4','1000';'EX_arab_D(e)','-10e-4','1000';'EX_cgly(e)','-10e-4','1000';'EX_cit(e)','-10e-4','1000';'EX_no2(e)','-10e-4','1000';'EX_no3(e)','-10e-4','1000';'EX_ocdca(e)','-10e-4','1000';'EX_orn(e)','-10e-4','1000';'EX_ptrc(e)','-10e-4','1000';'EX_sheme(e)','-10e-4','1000';'EX_spmd(e)','-10e-4','1000';'EX_mqn7(e)','-10e-4','1000';'EX_pydxn(e)','-10e-4','1000'};
western_diet={'EX_fru(e)','-0.14899','1000';'EX_glc_D(e)','-0.14899','1000';'EX_gal(e)','-0.14899','1000';'EX_man(e)','-0.14899','1000';'EX_mnl(e)','-0.14899','1000';'EX_fuc_L(e)','-0.14899','1000';'EX_glcn(e)','-0.14899','1000';'EX_rmn(e)','-0.14899','1000';'EX_arab_L(e)','-0.17878','1000';'EX_drib(e)','-0.17878','1000';'EX_rib_D(e)','-0.17878','1000';'EX_xyl_D(e)','-0.17878','1000';'EX_oxa(e)','-0.44696','1000';'EX_lcts(e)','-0.074493','1000';'EX_malt(e)','-0.074493','1000';'EX_sucr(e)','-0.074493','1000';'EX_melib(e)','-0.074493','1000';'EX_cellb(e)','-0.074493','1000';'EX_tre(e)','-0.074493','1000';'EX_strch1(e)','-0.25734','1000';'EX_amylopect900(e)','-1.5673e-05','1000';'EX_amylose300(e)','-4.7019e-05','1000';'EX_arabinan101(e)','-0.00016628','1000';'EX_arabinogal(e)','-2.1915e-05','1000';'EX_arabinoxyl(e)','-0.00030665','1000';'EX_bglc(e)','-7.05e-08','1000';'EX_cellul(e)','-2.8212e-05','1000';'EX_dextran40(e)','-0.00017632','1000';'EX_galmannan(e)','-1.4106e-05','1000';'EX_glcmannan(e)','-3.2881e-05','1000';'EX_homogal(e)','-0.00012823','1000';'EX_inulin(e)','-0.00047019','1000';'EX_kestopt(e)','-0.0028212','1000';'EX_levan1000(e)','-1.4106e-05','1000';'EX_lmn30(e)','-0.00047019','1000';'EX_lichn(e)','-8.2976e-05','1000';'EX_pect(e)','-3.3387e-05','1000';'EX_pullulan1200(e)','-1.1755e-05','1000';'EX_raffin(e)','-0.0047019','1000';'EX_rhamnogalurI(e)','-1.4492e-05','1000';'EX_rhamnogalurII(e)','-0.00026699','1000';'EX_starch1200(e)','-1.1755e-05','1000';'EX_xylan(e)','-3.2059e-05','1000';'EX_xyluglc(e)','-1.3146e-05','1000';'EX_arachd(e)','-0.003328','1000';'EX_chsterol(e)','-0.004958','1000';'EX_glyc(e)','-1.7997','1000';'EX_hdca(e)','-0.39637','1000';'EX_hdcea(e)','-0.036517','1000';'EX_lnlc(e)','-0.35911','1000';'EX_lnlnca(e)','-0.017565','1000';'EX_lnlncg(e)','-0.017565','1000';'EX_ocdca(e)','-0.16928','1000';'EX_ocdcea(e)','-0.68144','1000';'EX_octa(e)','-0.012943','1000';'EX_ttdca(e)','-0.068676','1000';'EX_ala_L(e)','-1','1000';'EX_cys_L(e)','-1','1000';'EX_ser_L(e)','-1','1000';'EX_arg_L(e)','-0.15','1000';'EX_his_L(e)','-0.15','1000';'EX_ile_L(e)','-0.15','1000';'EX_leu_L(e)','-0.15','1000';'EX_lys_L(e)','-0.15','1000';'EX_asn_L(e)','-0.225','1000';'EX_asp_L(e)','-0.225','1000';'EX_thr_L(e)','-0.225','1000';'EX_glu_L(e)','-0.18','1000';'EX_met_L(e)','-0.18','1000';'EX_gln_L(e)','-0.18','1000';'EX_pro_L(e)','-0.18','1000';'EX_val_L(e)','-0.18','1000';'EX_phe_L(e)','-1','1000';'EX_tyr_L(e)','-1','1000';'EX_gly(e)','-0.45','1000';'EX_trp_L(e)','-0.08182','1000';'EX_12dgr180(e)','-1','1000';'EX_26dap_M(e)','-1','1000';'EX_2dmmq8(e)','-1','1000';'EX_2obut(e)','-1','1000';'EX_3mop(e)','-1','1000';'EX_4abz(e)','-1','1000';'EX_4hbz(e)','-1','1000';'EX_5aop(e)','-1','1000';'EX_ac(e)','-1','1000';'EX_acald(e)','-1','1000';'EX_acgam(e)','-1','1000';'EX_acmana(e)','-1','1000';'EX_acnam(e)','-1','1000';'EX_ade(e)','-1','1000';'EX_adn(e)','-1','1000';'EX_adocbl(e)','-1','1000';'EX_akg(e)','-1','1000';'EX_ala_D(e)','-1','1000';'EX_amet(e)','-1','1000';'EX_amp(e)','-1','1000';'EX_anth(e)','-1','1000';'EX_arab_D(e)','-1','1000';'EX_avite1(e)','-1','1000';'EX_btn(e)','-1','1000';'EX_ca2(e)','-1','1000';'EX_cbl1(e)','-1','1000';'EX_cgly(e)','-1','1000';'EX_chol(e)','-1','1000';'EX_chor(e)','-1','1000';'EX_cit(e)','-1','1000';'EX_cl(e)','-1','1000';'EX_cobalt2(e)','-1','1000';'EX_csn(e)','-1','1000';'EX_cu2(e)','-1','1000';'EX_cytd(e)','-1','1000';'EX_dad_2(e)','-1','1000';'EX_dcyt(e)','-1','1000';'EX_ddca(e)','-1','1000';'EX_dgsn(e)','-1','1000';'EX_etoh(e)','-1','1000';'EX_fald(e)','-1','1000';'EX_fe2(e)','-1','1000';'EX_fe3(e)','-1','1000';'EX_fe3dcit(e)','-1','1000';'EX_fol(e)','-1','1000';'EX_for(e)','-1','1000';'EX_fum(e)','-1','1000';'EX_gam(e)','-1','1000';'EX_glu_D(e)','-1','1000';'EX_glyc3p(e)','-1','1000';'EX_gsn(e)','-1','1000';'EX_gthox(e)','-1','1000';'EX_gthrd(e)','-1','1000';'EX_gua(e)','-1','1000';'EX_h(e)','-1','1000';'EX_h2(e)','-1','1000';'EX_h2s(e)','-1','1000';'EX_hom_L(e)','-1','1000';'EX_hxan(e)','-1','1000';'EX_indole(e)','-1','1000';'EX_ins(e)','-1','1000';'EX_k(e)','-1','1000';'EX_lac_L(e)','-1','1000';'EX_lanost(e)','-1','1000';'EX_mal_L(e)','-1','1000';'EX_metsox_S_L(e)','-1','1000';'EX_mg2(e)','-1','1000';'EX_mn2(e)','-1','1000';'EX_mobd(e)','-1','1000';'EX_mqn7(e)','-1','1000';'EX_mqn8(e)','-1','1000';'EX_na1(e)','-1','1000';'EX_nac(e)','-1','1000';'EX_ncam(e)','-1','1000';'EX_nmn(e)','-1','1000';'EX_no2(e)','-1','1000';'EX_no2(e)','-1','1000';'EX_no3(e)','-1','1000';'EX_orn(e)','-1','1000';'EX_pheme(e)','-1','1000';'EX_phyQ(e)','-1','1000';'EX_pi(e)','-1','1000';'EX_pime(e)','-1','1000';'EX_pnto_R(e)','-1','1000';'EX_ptrc(e)','-1','1000';'EX_pydam(e)','-1','1000';'EX_pydx(e)','-1','1000';'EX_pydx5p(e)','-1','1000';'EX_pydxn(e)','-1','1000';'EX_q8(e)','-1','1000';'EX_retinol(e)','-1','1000';'EX_ribflv(e)','-1','1000';'EX_sel(e)','-1','1000';'EX_sheme(e)','-1','1000';'EX_so4(e)','-1','1000';'EX_spmd(e)','-1','1000';'EX_succ(e)','-1','1000';'EX_thf(e)','-1','1000';'EX_thm(e)','-1','1000';'EX_thymd(e)','-1','1000';'EX_ura(e)','-1','1000';'EX_uri(e)','-1','1000';'EX_vitd3(e)','-1','1000';'EX_xan(e)','-1','1000';'EX_zn2(e)','-1','1000';'EX_meoh(e)','-10','1000';'EX_h2o(e)','-10','1000'};
ycfa_media={'EX_pi(e)','-7.605312000000001','1000';'EX_k(e)','-15.22291','1000';'EX_4abz(e)','-0.00021899999999999998','1000';'EX_ac(e)','-33.0','1000';'EX_arg_L(e)','-1.9804739999999998','1000';'EX_arg_D(e)','-1.9804739999999998','1000';'EX_btn(e)','-5.47e-05','1000';'EX_ca(e)','-0.8109569999999999','1000';'EX_cl(e)','-24.457700999999997','1000';'EX_cbl1(e)','-7.38e-06','1000';'EX_cys_L(e)','-6.716115000000001','1000';'EX_h(e)','-6.344423','1000';'EX_so4(e)','-3.8777579999999996','1000';'EX_h2o(e)','-27754.6489','1000';'EX_fol(e)','-0.00011300000000000001','1000';'EX_na1(e)','-78.695299','1000';'EX_glc_D(e)','-27.753107999999997','1000';'EX_mg2(e)','-2.243342','1000';'EX_mn2(e)','-0.01','1000';'EX_hco3(e)','-47.613378999999995','1000';'EX_nh4(e)','-1.8318919999999999','1000';'EX_fe3(e)','-0.0218','1000';'EX_pydam(e)','-0.000622','1000';'EX_thm(e)','-0.00104','1000';'EX_ppa(e)','-9.0','1000';'EX_isobut(e)','-1.0','1000';'EX_isoval(e)','-1.0','1000';'EX_M03134(e)','-1.0','1000';'EX_ala_D(e)','-4.377525','1000';'EX_ala_L(e)','-4.377525','1000';'EX_asp_L(e)','-5.327528','1000';'EX_gln_L(e)','-0.17105499999999998','1000';'EX_gln_D(e)','-0.17105499999999998','1000';'EX_glu_D(e)','-6.353279','1000';'EX_glu_L(e)','-6.353279','1000';'EX_gly(e)','-3.796496','1000';'EX_his_L(e)','-0.918456','1000';'EX_ile_L(e)','-2.706389','1000';'EX_leu_L(e)','-4.059588000000001','1000';'EX_lys_L(e)','-2.770368','1000';'EX_met_L(e)','-0.469139','1000';'EX_phe_L(e)','-2.270109','1000';'EX_pro_L(e)','-2.1714510000000002','1000';'EX_ser_L(e)','-1.807963','1000';'EX_asn_L(e)','-1.112662','1000';'EX_trp_L(e)','-0.1591','1000';'EX_tyr_L(e)','-0.88305','1000';'EX_val_L(e)','-3.0516560000000004','1000';'EX_ca2(e)','-0.145336','1000';'EX_fe2(e)','-0.0218','1000';'EX_thr_L(e)','-0.33579600000000004','1000';'EX_cd2(e)','-3.34e-05','1000';'EX_cobalt2(e)','-0.01','1000';'EX_cu2(e)','-0.01','1000';'EX_ni2(e)','-0.000272','1000';'EX_zn2(e)','-0.141403','1000';'EX_ade(e)','-0.163546','1000';'EX_gua(e)','-0.15301499999999998','1000';'EX_csn(e)','-0.09809999999999999','1000';'EX_ura(e)','-0.140738','1000';'EX_man(e)','-0.63556','1000';'EX_pnto_R(e)','-0.0012100000000000001','1000';'EX_ascb_L(e)','-4.26e-05','1000';'EX_nac(e)','-0.00816','1000';'EX_pydx(e)','-0.000222','1000';'EX_chol(e)','-0.007679999999999999','1000';'EX_adocbl(e)','-1.2899999999999999e-09','1000';'EX_mobd(e)','-0.000154','1000';'EX_n2(e)','-800.0','1000';'EX_h2(e)','-100.0','1000';'EX_26dap_M(e)','-10e-4','1000';'EX_h2s(e)','-10e-4','1000';'EX_ribflv(e)','-10e-4','1000';'EX_thymd(e)','-10e-4','1000';'EX_12dgr180(e)','-10e-4','1000';'EX_hxan(e)','-10e-4','1000';'EX_ddca(e)','-10e-4','1000';'EX_ttdca(e)','-10e-4','1000';'EX_pheme(e)','-10e-4','1000';'EX_q8(e)','-10e-4','1000';'EX_arab_D(e)','-10e-4','1000';'EX_cgly(e)','-10e-4','1000';'EX_cit(e)','-10e-4','1000';'EX_no2(e)','-10e-4','1000';'EX_no3(e)','-10e-4','1000';'EX_ocdca(e)','-10e-4','1000';'EX_orn(e)','-10e-4','1000';'EX_ptrc(e)','-10e-4','1000';'EX_sheme(e)','-10e-4','1000';'EX_spmd(e)','-10e-4','1000';'EX_mqn7(e)','-10e-4','1000'};
LB_media={'EX_na1(e)', '-1.87E+02', '1000'; 'EX_cl(e)', '-1.72E+02', '1000'; 'EX_ca2(e)', '-8.01E-02', '1000'; 'EX_fe2(e)', '-4.56E-02', '1000'; 'EX_fe3(e)', '4.56E-02', '1000'; 'EX_so4(e)', '-3.90E-01', '1000'; 'EX_pi(e)', '-4.44E+00', '1000'; 'EX_mg2(e)', '-1.75E+00', '1000'; 'EX_k(e)', '-2.58E+00', '1000'; 'EX_ala_L(e)', '-6.73E+00', '1000'; 'EX_asp_L(e)', '-5.90E+00', '1000'; 'EX_asn_L(e)', '-8.33E-01', '1000'; 'EX_glu_L(e)', '-1.35E+01', '1000'; 'EX_gln_L(e)', '-1.37E-01', '1000'; 'EX_gly(e)', '-4.26E+00', '1000'; 'EX_his_L(e)', '-1.64E+00', '1000'; 'EX_ile_L(e)', '-5.34E+00', '1000'; 'EX_leu_L(e)', '-7.28E+00', '1000'; 'EX_lys_L(e)', '-5.81E+00', '1000'; 'EX_met_L(e)', '-1.68E+00', '1000'; 'EX_phe_L(e)', '-3.93E+00', '1000'; 'EX_pro_L(e)','-6.60E+00','1000';'EX_ser_L(e)','-2.85E+00','1000';'EX_thr_L(e)','-2.18E+00','1000';'EX_trp_L(e)','-5.14E-01','1000';'EX_tyr_L(e)','-1.05E+00','1000';'EX_val_L(e)','-6.53E+00','1000';'EX_arg_L(e)','-3.62E+00','1000';'EX_cys_L(e)','-3.33E-01','1000';'EX_cd2(e)','-6.67E-05','1000';'EX_cobalt2(e)','-2.97E-04','1000';'EX_cu(e)','-5.29E-03','1000';'EX_cu2(e)','-5.29E-03','1000';'EX_mn2(e)','-1.10E-02','1000';'EX_ni2(e)','-2.08E-03','1000';'EX_zn2(e)','-4.94E-01','1000';'EX_ade(e)','-3.27E-01','1000';'EX_gua(e)','-3.06E-01','1000';'EX_csn(e)','1.96E-01','1000';'EX_ura(e)','-2.81E-01','1000';'EX_nh4(e)','-2.30E+00','1000';'EX_man(e)','-1.27E+00','1000';'EX_pnto_R(e)','-2.42E-03','1000';'EX_btn(e)','-2.76E-05','1000';'EX_ascb_L(e)','-8.52E-05','1000';'EX_thm(e)','-2.07E-03','1000';'EX_nac(e)','-1.63E-02','1000';'EX_pydx(e)','-4.43E-04','1000';'EX_chol(e)','-1.54E-02','1000';'EX_adocbl(e)','-2.58E-09','1000';'EX_o2(e)','-1.82E+01','1000';'EX_h2o(e)','-5.55E+04','1000';'EX_h(e)','-1.00E-04','1000';'EX_mobd(e)','-3.07E-04','1000'};


pyg_altered = {'EX_n2(e)','-100','1000';'EX_na1(e)','-100','1000';'EX_cl(e)','-100','1000';'EX_ca2(e)','-100','1000';'EX_fe2(e)','-100','1000';'EX_mg2(e)','-100','1000';'EX_k(e)','-100','1000';'EX_so4(e)','-100','1000';'EX_pi(e)','-100','1000';'EX_ala_L(e)','-100','1000';'EX_arg_L(e)','-100','1000';'EX_asn_L(e)','-100','1000';'EX_asp_L(e)','-100','1000';'EX_Lcystin(e)','-100','1000';'EX_glu_L(e)','-100','1000';'EX_gln_L(e)','-100','1000';'EX_gly(e)','-100','1000';'EX_his_L(e)','-100','1000';'EX_ile_L(e)','-100','1000';'EX_leu_L(e)','-100','1000';'EX_lys_l(e)','-100','1000';'EX_met_L(e)','-100','1000';'EX_phe_L(e)','-100','1000';'EX_pro_L(e)','-100','1000';'EX_ser_L(e)','-100','1000';'EX_thr_L(e)','-100','1000';'EX_trp_L(e)','-100','1000';'EX_tyr_L(e)','-100','1000';'EX_val_L(e)','-100','1000';'EX_glc_D(e)','-100','1000';'EX_hco3(e)','-100','1000';'EX_phyQ(e)','-100','1000';'EX_etoh(e)','-100','1000';'EX_pheme(e)','-100','1000';'EX_oh1(e)','-100','1000';'EX_cys_L(e)','-100','1000';'EX_M02144(e)','-100','1000';'EX_h2o(e)','-55013.2623','1000';'EX_h(e)','-6.30957E-05','1000'};

%%

p_copri_alt_pyg_alt = useDiet(p_copri_model_alt, pyg_altered);

p_copri_alt_pyg_alt= changeRxnBounds(p_copri_alt_pyg_alt, 'biomass525', 0.001, 'l');

test_sol_pyg_alt = optimizeCbModel(p_copri_alt_pyg_alt, 'max');


%%

% Initialize a logical vector to store the presence of "EX_" reactions
isExchangeReaction = false(length(p_copri_model_alt.rxns), 1);

% Loop through the reactions in the model
for i = 1:length(p_copri_model_alt.rxns)
    % Check if the reaction ID is present in pyg_media
    if any(strcmp(pyg_media(:, 1), p_copri_model_alt.rxns{i}))
        isExchangeReaction(i) = true;
    else
        isExchangeReaction(i) = false;
    end
end

% Display the vector of exchange reaction presence
disp('Vector indicating the presence of "EX_" reactions in pyg_media:');
disp(isExchangeReaction);

excludedReactionLB = find(isExchangeReaction); % indices of exchange reactions
excludedReactionLB = p_copri_model_alt.lb(excludedReactionLB); % lower bounds


%% Set diet/media for models

% use diet function not working for some reason 


f_prau_pyg = useDiet(f_prau_model, pyg_test);
f_prau_ABB = useDiet(f_prau_model, ABB_media);


p_copri_alt_LB = useDiet(p_copri_model_alt, LB_media);
f_prau_LB = useDiet(f_prau_model, LB_media);

p_copri_pyg = useDiet(p_copri_model, pyg_test);
p_copri_alt_pyg = useDiet(p_copri_model_alt, pyg_test);
combined_model_pyg = useDiet(combined_model, pyg_media);
combined_model_abb = useDiet(combined_model, ABB_media);

p_copri_alt_western = useDiet(p_copri_model_alt, western_diet);


%% Just trying to simulate growth of P. copri on pyg 

p_copri_alt_pyg = useDiet(p_copri_model_alt, pyg_test);
p_copri_alt_abb = useDiet(p_copri_model_alt, ABB_media);
p_copri_alt_ycfa = useDiet(p_copri_model_alt, ycfa_media);


p_copri_alt_pyg = changeRxnBounds(p_copri_alt_pyg, 'biomass525', 0.001, 'l');
p_copri_alt_abb = changeRxnBounds(p_copri_alt_abb, 'biomass525', 0.001, 'l');
p_copri_alt_ycfa = changeRxnBounds(p_copri_alt_ycfa, 'biomass525', 0.001, 'l');
p_copri_alt_LB= changeRxnBounds(p_copri_alt_LB, 'biomass525', 0.001, 'l');

test_sol = optimizeCbModel(p_copri_alt_pyg, 'max');
test_abb = optimizeCbModel(p_copri_alt_abb);
test_ycga = optimizeCbModel(p_copri_alt_ycga);


%%

param = struct();
param.internalRelax = 2;  % Allow to relax bounds on all internal reactions
param.exchangeRelax = 0;  % Do not allow to relax bounds on exchange reactions

% Define the list of excluded reactions for lower bound relaxation
excludedReactionIDs = {'biomass525'};

% Create a logical vector to exclude these reactions from relaxation
excludedReactionLB = false(size(p_copri_alt_pyg.rxns)); % Initialize as false
for i = 1:length(p_copri_alt_pyg.rxns)
    % Check if the current reaction is in the excludedReactionIDs list
    if any(strcmp(p_copri_alt_pyg.rxns{i}, excludedReactionIDs))
        excludedReactionLB(i) = true; % Exclude this reaction from relaxation
    end
end
    % Find the index of the reaction in the model
    %rxnIndex = find(strcmp(p_copri_alt_pyg.rxns, excludedReactionIDs{i}));
    %if ~isempty(rxnIndex)z
    %    excludedReactionLB(rxnIndex) = true; % Set the corresponding index to true
    %end


%%
% Set all lower bounds to zero
p_copri_alt_pyg.lb(:) = 0;


% Iterate through the dietary constraints and set bounds directly
for i = 1:size(pyg_test, 1)
    % Extract reaction ID, lower bound, and upper bound
    rxnID = pyg_test{i, 1};
    lb = str2double(pyg_test{i, 2});
    ub = str2double(pyg_test{i, 3});
    
    % Set the bounds for the reaction in the model
    p_copri_alt_pyg = changeRxnBounds(p_copri_alt_pyg, rxnID, lb, 'l'); % Lower bound
    p_copri_alt_pyg = changeRxnBounds(p_copri_alt_pyg, rxnID, ub, 'u'); % Upper bound
end




%%
%p_copri_alt_pyg = useDiet(p_copri_model_alt, pyg_test);
p_copri_alt_pyg = changeRxnBounds(p_copri_alt_pyg, 'biomass525', 0.001, 'l');

%%

% Add the excludedReactionLB to the param structure
param.excludedReactionLB = excludedReactionLB;


test_relax = relaxedFBA(p_copri_alt_pyg, param);

%%


p_copri_alt_pyg.lb = p_copri_alt_pyg.lb - test_relax.p;

%%


p_copri_alt_pyg_sol = optimizeCbModel(p_copri_alt_pyg);


%% set biomass as objective function and set minimal growth

% for each paired model, set both biomass objective functions as
% objectives
%biomass1 = strcat(info{i, 2}, '_', info{i, 3});
%biomass2 = strcat(info{i, 4}, '_', info{i, 5});
%model1biomass = find(ismember(pairedModel.rxns, biomass1));
%pairedModel.c(model1biomass, 1) = 1;
%model2biomass = find(ismember(pairedModel.rxns, biomass2));
%pairedModel.c(model2biomass, 1) = 1;
% enforce minimal growth
f_prau_model = changeRxnBounds(f_prau_model, 'bio1', 0.001, 'l');

f_prau_pyg = changeRxnBounds(f_prau_pyg, 'bio1', 0.001, 'l');
p_copri_pyg = changeRxnBounds(p_copri_pyg, 'bio1', 0.001, 'l');
p_copri_alt_pyg = changeRxnBounds(p_copri_alt_pyg, 'biomass525', 0.001, 'l');
p_copri_alt_western = changeRxnBounds(p_copri_alt_western, 'biomass525', 0.001, 'l');

f_prau_LB = changeRxnBounds(f_prau_LB, 'bio1', 0.001, 'l');

combined_model_pyg.c('Prevotella_copri_ERR1022397_bio1', 1) = 1;
combined_model_pyg.c('Faecalibacterium_prausnitzii_ERR1022484_bio1', 1) = 1;
% enforce minimal growth
combined_model_pyg = changeRxnBounds(combined_model_abb, 'Prevotella_copri_ERR1022397_bio1', 0.001, 'l');
combined_model_pyg = changeRxnBounds(combined_model_abb, 'Faecalibacterium_prausnitzii_ERR1022484_bio1', 0.001, 'l');

combined_model_abb.c('Prevotella_copri_ERR1022397_bio1', 1) = 1;
combined_model_abb.c('Faecalibacterium_prausnitzii_ERR1022484_bio1', 1) = 1;
% enforce minimal growth
combined_model_abb = changeRxnBounds(combined_model_abb, 'Prevotella_copri_ERR1022397_bio1', 0.001, 'l');
combined_model_abb = changeRxnBounds(combined_model_abb, 'Faecalibacterium_prausnitzii_ERR1022484_bio1', 0.001, 'l');





%% Optimize mono culture growth 

f_prau_pyg_sol = optimizeCbModel(f_prau_pyg);
p_copri_pyg_sol = optimizeCbModel(p_copri_pyg);

p_copri_alt_pyg_sol = optimizeCbModel(p_copri_alt_pyg);

combined_pyg_sol = optimizeCbModel(combined_model_pyg);
combined_abb_sol = optimizeCbModel(combined_model_abb);

f_prau_sol = optimizeCbModel(f_prau_model);

p_copri_alt_west_sol = optimizeCbModel(p_copri_alt_western);




%% Examine solution

LB_test = optimizeCbModel(p_copri_alt_LB);

f_prau_LB_sol = optimizeCbModel(f_prau_LB)


