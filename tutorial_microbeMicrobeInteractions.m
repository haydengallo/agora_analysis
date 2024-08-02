%% Computation and analysis of microbe-microbe metabolic interactions
%% Author: Almut Heinken, Molecular Systems Physiology Group, National University of Ireland Galway.
% This tutorial demonstrates how to join a given list of microbial COBRA models 
% in all possible combinations and compute the metabolic
% 
% interactions between the microbes depending on the implemented diet. Moreover, 
% the tradeoff between the growth of different joined microbes is computed. The 
% tutorial can be adapted to any number of AGORA models and dietary conditions 
% analyzed.
%% Requirements
% This tutorial requires the Parallel Computing Toolbox, Bioinformatics Toolbox, 
% and Statistics and Machine Learning Toolbox add-ons in MATLAB.
%% Initialize the COBRA Toolbox

initCobraToolbox
%% Prepare input data and models
% We will use the AGORA resource (Magnusdottir et al., Nat Biotechnol. 2017 
% Jan;35(1):81-89) in this tutorial. AGORA version 1.03 is available at https://github.com/VirtualMetabolicHuman/AGORA. 
% Download AGORA and place the models into a folder.

%websave('AGORA-master.zip','https://github.com/VirtualMetabolicHuman/AGORA/archive/master.zip')
%try
%    unzip('AGORA-master')
%end
%modPath = [pwd filesep 'AGORA-master' filesep 'CurrentVersion' filesep 'AGORA_1_03' filesep' 'AGORA_1_03_mat'];
modPath = '/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/';
%% 
% Import a file with information on the AGORA organisms including reconstruction 
% names and taxonomy.

[~,infoFile,~]=xlsread('/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/p_copri_f_prau_only.xlsx');
%[~,infoFile,~]=xlsread('/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/AGORA2_infoFile_full.xlsx');
%% 
% Note: if you get a 'file not found' error for AGORA_infoFile.xlsx, please 
% run initCobraToolbox and ensure that the cobratoolbox/papers folder is in the 
% MATLAB path.
%% Creation of pairwise models
% For the sake of this tutorial, we will use ten random AGORA reconstructions 
% from the info file.

%modelList = infoFile(randi([2 length(infoFile)],1,10),1);
%% 
% Uncomment the following line to join all AGORA reconstructions in all combinations. 
% NOTE: this is very time-consuming due to the large number of model combinations 
% analyzed.
% 
modelList=infoFile(2:end,1);
% 
% You may also enter a custom selection of AGORA reconstructions as a cell array 
% named modelList.
% 
% Let us define some parameters for joining the models. Set the coupling factor 
% c, which defined how the flux through all reactions in a model is coupled to 
% the flux through its biomass reaction. Allowed flux span through each reaction= 
% -(c * flux(biomass)) to +(c * flux(biomass)).

c = 400;
%% 
% Set the threshold u, which defines the flux through each reaction that is 
% allowed if flux through the biomass reaction is zero.

u = 0;
%% 
% Define whether or not genes from the models are merged and kept in the joined 
% models. If set to true the joining is more time-consuming.

mergeGenes = false;
%% 
% Define the number workers for parallel pool to allow parallel computing.Recommended 
% if a large number of microbe models is computed. Set to zero if parallel computing 
% is not available.

numWorkers = 10;
%% 
% path where to save results

mkdir('Results_p_copri_f_prau');
resPath = [pwd filesep 'Results_p_copri_f_prau'];

%mkdir('interactions_test');
%resPath = [pwd filesep 'interactions_test'];

%% 
% Join the models in all possible combinations.

joinModelsPairwiseFromList(modelList,modPath,'pairwiseModelFolder', resPath,'c',c,'u',u,'mergeGenesFlag',mergeGenes,'numWorkers',numWorkers);
%% Computation of pairwise interactions
% The interactions between all microbes joined in the first step will be simulated 
% on given dietary conditions. Here, we will use four dietary conditions used 
% in Magnusdottir et al., Nat Biotechnol. 2017.: Western Diet without oxygen, 
% Western Diet with oxygen, High fiber diet without oxygen, and High fiber diet 
% with oxygen. Let us define the input parameters for the simulation of pairwise 
% interactions.
% 
% Name the four dietary conditions that will be simulated.

%conditions = {'PYG','test_ABB','test_YCAG','test_YCGMS','test_YCFA','test_mMCB'};

conditions = {'western','PYG_with_all'}%'test_ABB','test_YCAG','LB'}%,'test_YCGMS','test_YCFA','test_mMCB'};



%% 
% Define the corresponding constraints to implement for each diet. The input 
% file needs to be a string array.



% ABB

%dietConstraints{2}={'EX_pi[u]','-3.777186','1000';'EX_k[u]','-12.084479','1000';'EX_ac[u]','-0.0','1000';'EX_arg_L[u]','-9.357036999999998','1000';'EX_arg_D[u]','-9.357036999999998','1000';'EX_btn[u]','-3.86e-05','1000';'EX_cl[u]','-90.79522399999998','1000';'EX_cbl1[u]','-10e-4','1000';'EX_co2[u]','-100.0','1000';'EX_cys_L[u]','-3.8169380000000004','1000';'EX_h[u]','-3.172326','1000';'EX_dtt[u]','-6.482982','1000';'EX_so4[u]','-3.9465080000000006','1000';'EX_h2o[u]','-77713.01693','1000';'EX_na1[u]','-131.369616','1000';'EX_glc_D[u]','-5.550622','1000';'EX_h2[u]','-50.0','1000';'EX_mg2[u]','-3.393075','1000';'EX_mqn8[u]','-0.0011200000000000001','1000';'EX_mn2[u]','-0.01','1000';'EX_n2[u]','-850.0','1000';'EX_hco3[u]','-4.761338','1000';'EX_pyr[u]','-9.087605','1000';'EX_succ[u]','-1.850892','1000';'EX_nh4[u]','-3.2222280000000003','1000';'EX_ppi[u]','-0.6709520000000001','1000';'EX_fe3[u]','-0.731952','1000';'EX_strch1[u]','-0.504602','1000';'EX_thiog[u]','-4.382121','1000';'EX_thm[u]','-0.0029000000000000002','1000';'EX_ala_D[u]','-8.889688999999999','1000';'EX_ala_L[u]','-8.889688999999999','1000';'EX_asp_L[u]','-9.718608','1000';'EX_gln_L[u]','-0.314768','1000';'EX_gln_D[u]','-0.314768','1000';'EX_glu_D[u]','-12.081918','1000';'EX_glu_L[u]','-12.081918','1000';'EX_gly[u]','-7.273324000000001','1000';'EX_his_L[u]','-1.7208890000000001','1000';'EX_ile_L[u]','-5.01633','1000';'EX_leu_L[u]','-7.433019','1000';'EX_lys_L[u]','-5.376564999999999','1000';'EX_met_L[u]','-0.911472','1000';'EX_phe_L[u]','-4.104351','1000';'EX_pro_L[u]','-3.9954629999999995','1000';'EX_ser_L[u]','-3.3494800000000002','1000';'EX_asn_L[u]','-2.0073269999999996','1000';'EX_trp_L[u]','-0.328062','1000';'EX_tyr_L[u]','-1.6115650000000001','1000';'EX_val_L[u]','-5.778952','1000';'EX_ca2[u]','-0.242261','1000';'EX_fe2[u]','-0.061','1000';'EX_thr_L[u]','-0.940228','1000';'EX_cd2[u]','-9.340000000000001e-05','1000';'EX_cobalt2[u]','-0.01','1000';'EX_cu2[u]','-0.01','1000';'EX_ni2[u]','-0.000763','1000';'EX_zn2[u]','-0.395929','1000';'EX_ade[u]','-0.457929','1000';'EX_gua[u]','-0.428442','1000';'EX_csn[u]','-0.274698','1000';'EX_ura[u]','-0.394065','1000';'EX_man[u]','-1.779569','1000';'EX_pnto_R[u]','-0.0033799999999999998','1000';'EX_ascb_L[u]','-0.00011899999999999999','1000';'EX_nac[u]','-0.0229','1000';'EX_pydx[u]','-0.0006209999999999999','1000';'EX_chol[u]','-0.021500000000000002','1000';'EX_adocbl[u]','-3.62e-09','1000';'EX_mobd[u]','-0.00043','1000';'EX_26dap_M[u]','-10e-4','1000';'EX_h2s[u]','-10e-4','1000';'EX_ribflv[u]','-10e-4','1000';'EX_thymd[u]','-10e-4','1000';'EX_fol[u]','-10e-4','1000';'EX_12dgr180[u]','-10e-4','1000';'EX_hxan[u]','-10e-4','1000';'EX_ddca[u]','-10e-4','1000';'EX_ttdca[u]','-10e-4','1000';'EX_pheme[u]','-10e-4','1000';'EX_q8[u]','-10e-4','1000';'EX_arab_D[u]','-10e-4','1000';'EX_cgly[u]','-10e-4','1000';'EX_cit[u]','-10e-4','1000';'EX_no2[u]','-10e-4','1000';'EX_no3[u]','-10e-4','1000';'EX_ocdca[u]','-10e-4','1000';'EX_orn[u]','-10e-4','1000';'EX_ptrc[u]','-10e-4','1000';'EX_sheme[u]','-10e-4','1000';'EX_spmd[u]','-10e-4','1000';'EX_mqn7[u]','-10e-4','1000';'EX_pydxn[u]','-10e-4','1000'};

% pyg with metals
%dietConstraints{2} = {'EX_zn2[u]','-1','1000';'EX_mn2[u]','-1','1000';'EX_fe3[u]','-1','1000';'EX_cu2[u]','-1','1000';'EX_cobalt2[u]','-1','1000';'EX_n2[u]','-227.8571429','1000';'EX_na1[u]','-20.51163798','1000';'EX_cl[u]','-5.941065976','1000';'EX_ca2[u]','-0.173941043','1000';'EX_fe2[u]','-0.016053362','1000';'EX_mg2[u]','-0.474477191','1000';'EX_k[u]','-35.39748582','1000';'EX_so4[u]','-0.710983412','1000';'EX_pi[u]','-18.29826648','1000';'EX_ala_L[u]','-16.89227108','1000';'EX_arg_L[u]','-5.338568575','1000';'EX_asn_L[u]','-1.286718791','1000';'EX_asp_L[u]','-10.81893313','1000';'EX_Lcystin[u]','-0.187272152','1000';'EX_glu_L[u]','-18.56754922','1000';'EX_gln_L[u]','-0.205274178','1000';'EX_gly[u]','-20.3151851','1000';'EX_his_L[u]','-3.319218598','1000';'EX_ile_L[u]','-8.195159139','1000';'EX_leu_L[u]','-11.2826377','1000';'EX_lys_l[u]','-9.884397018','1000';'EX_met_L[u]','-2.144657123','1000';'EX_phe_L[u]','-6.083829725','1000';'EX_pro_L[u]','-11.89938505','1000';'EX_ser_L[u]','-4.424652451','1000';'EX_thr_L[u]','-3.567830759','1000';'EX_trp_L[u]','-0.685504997','1000';'EX_tyr_L[u]','-1.683306566','1000';'EX_val_L[u]','-10.37149589','1000';'EX_glc_D[u]','-27.7537245498346','1000';'EX_hco3[u]','-1.190379826','1000';'EX_phyQ[u]','-0.002172408','1000';'EX_etoh[u]','-3.237913066','1000';'EX_pheme[u]','-0.0076693','1000';'EX_oh1[u]','-0.099987502','1000';'EX_cys_L[u]','-2.846894039','1000';'EX_M02144[u]','-2.846894039','1000';'EX_h2o[u]','-55013.2623','1000';'EX_h[u]','-6.30957E-05','1000'};

% pyg with addition for both P. copri and F. prau
dietConstraints{2} = {'EX_q8[u]','-1','1000';'EX_pydx[u]','-1','1000';'EX_nac[u]','-1','1000';'EX_lys_L[u]','-1','1000';'EX_hxan[u]','-1','1000';'EX_ade[u]','-1','1000';'EX_thymd[u]','-1','1000';'EX_thm[u]','-1','1000';'EX_ribflv[u]','-1','1000';'EX_pnto_R[u]','-1','1000';'EX_nac[u]','-1','1000';'EX_fol[u]','-1','1000';'EX_zn2[u]','-1','1000';'EX_mn2[u]','-1','1000';'EX_fe3[u]','-1','1000';'EX_cu2[u]','-1','1000';'EX_cobalt2[u]','-1','1000';'EX_n2[u]','-227.8571429','1000';'EX_na1[u]','-20.51163798','1000';'EX_cl[u]','-5.941065976','1000';'EX_ca2[u]','-0.173941043','1000';'EX_fe2[u]','-0.016053362','1000';'EX_mg2[u]','-0.474477191','1000';'EX_k[u]','-35.39748582','1000';'EX_so4[u]','-0.710983412','1000';'EX_pi[u]','-18.29826648','1000';'EX_ala_L[u]','-16.89227108','1000';'EX_arg_L[u]','-5.338568575','1000';'EX_asn_L[u]','-1.286718791','1000';'EX_asp_L[u]','-10.81893313','1000';'EX_Lcystin[u]','-0.187272152','1000';'EX_glu_L[u]','-18.56754922','1000';'EX_gln_L[u]','-0.205274178','1000';'EX_gly[u]','-20.3151851','1000';'EX_his_L[u]','-3.319218598','1000';'EX_ile_L[u]','-8.195159139','1000';'EX_leu_L[u]','-11.2826377','1000';'EX_lys_l[u]','-9.884397018','1000';'EX_met_L[u]','-2.144657123','1000';'EX_phe_L[u]','-6.083829725','1000';'EX_pro_L[u]','-11.89938505','1000';'EX_ser_L[u]','-4.424652451','1000';'EX_thr_L[u]','-3.567830759','1000';'EX_trp_L[u]','-0.685504997','1000';'EX_tyr_L[u]','-1.683306566','1000';'EX_val_L[u]','-10.37149589','1000';'EX_glc_D[u]','-27.7537245498346','1000';'EX_hco3[u]','-1.190379826','1000';'EX_phyQ[u]','-0.002172408','1000';'EX_etoh[u]','-3.237913066','1000';'EX_pheme[u]','-0.0076693','1000';'EX_oh1[u]','-0.099987502','1000';'EX_cys_L[u]','-2.846894039','1000';'EX_M02144[u]','-2.846894039','1000';'EX_h2o[u]','-55013.2623','1000';'EX_h[u]','-6.30957E-05','1000'};




%Other
%dietConstraints{2}={'EX_na1[u]', '-1.87E+02', '1000'; 'EX_cl[u]', '-1.72E+02', '1000'; 'EX_ca2[u]', '-8.01E-02', '1000'; 'EX_fe2[u]', '-4.56E-02', '1000'; 'EX_fe3[u]', '4.56E-02', '1000'; 'EX_so4[u]', '-3.90E-01', '1000'; 'EX_pi[u]', '-4.44E+00', '1000'; 'EX_mg2[u]', '-1.75E+00', '1000'; 'EX_k[u]', '-2.58E+00', '1000'; 'EX_ala_L[u]', '-6.73E+00', '1000'; 'EX_asp_L[u]', '-5.90E+00', '1000'; 'EX_asn_L[u]', '-8.33E-01', '1000'; 'EX_glu_L[u]', '-1.35E+01', '1000'; 'EX_gln_L[u]', '-1.37E-01', '1000'; 'EX_gly[u]', '-4.26E+00', '1000'; 'EX_his_L[u]', '-1.64E+00', '1000'; 'EX_ile_L[u]', '-5.34E+00', '1000'; 'EX_leu_L[u]', '-7.28E+00', '1000'; 'EX_lys_L[u]', '-5.81E+00', '1000'; 'EX_met_L[u]', '-1.68E+00', '1000'; 'EX_phe_L[u]', '-3.93E+00', '1000'; 'EX_pro_L[u]','-6.60E+00','1000';'EX_ser_L[u]','-2.85E+00','1000';'EX_thr_L[u]','-2.18E+00','1000';'EX_trp_L[u]','-5.14E-01','1000';'EX_tyr_L[u]','-1.05E+00','1000';'EX_val_L[u]','-6.53E+00','1000';'EX_arg_L[u]','-3.62E+00','1000';'EX_cys_L[u]','-3.33E-01','1000';'EX_cd2[u]','-6.67E-05','1000';'EX_cobalt2[u]','-2.97E-04','1000';'EX_cu[u]','-5.29E-03','1000';'EX_cu2[u]','-5.29E-03','1000';'EX_mn2[u]','-1.10E-02','1000';'EX_ni2[u]','-2.08E-03','1000';'EX_zn2[u]','-4.94E-01','1000';'EX_ade[u]','-3.27E-01','1000';'EX_gua[u]','-3.06E-01','1000';'EX_csn[u]','1.96E-01','1000';'EX_ura[u]','-2.81E-01','1000';'EX_nh4[u]','-2.30E+00','1000';'EX_man[u]','-1.27E+00','1000';'EX_pnto_R[u]','-2.42E-03','1000';'EX_btn[u]','-2.76E-05','1000';'EX_ascb_L[u]','-8.52E-05','1000';'EX_thm[u]','-2.07E-03','1000';'EX_nac[u]','-1.63E-02','1000';'EX_pydx[u]','-4.43E-04','1000';'EX_chol[u]','-1.54E-02','1000';'EX_adocbl[u]','-2.58E-09','1000';'EX_o2[u]','-1.82E+01','1000';'EX_h2o[u]','-5.55E+04','1000';'EX_h[u]','-1.00E-04','1000';'EX_mobd[u]','-3.07E-04','1000'};

%PYG
%dietConstraints{1}={'EX_mobd[u]','-3.07E-04','1000';'EX_ac[u]','-0.0','1000';'EX_n2[u]','-227.8571429','1000';'EX_na1[u]','-20.51163798','1000';'EX_cl[u]','-5.941065976','1000';'EX_ca2[u]','-0.173941043','1000';'EX_fe2[u]','-0.016053362','1000';'EX_mg2[u]','-0.474477191','1000';'EX_k[u]','-35.39748582','1000';'EX_so4[u]','-0.710983412','1000';'EX_pi[u]','-18.29826648','1000';'EX_ala_L[u]','-16.89227108','1000';'EX_arg_L[u]','-5.338568575','1000';'EX_asn_L[u]','-1.286718791','1000';'EX_asp_L[u]','-10.81893313','1000';'EX_Lcystin[u]','-0.187272152','1000';'EX_glu_L[u]','-18.56754922','1000';'EX_gln_L[u]','-0.205274178','1000';'EX_gly[u]','-20.3151851','1000';'EX_his_L[u]','-3.319218598','1000';'EX_ile_L[u]','-8.195159139','1000';'EX_leu_L[u]','-11.2826377','1000';'EX_lys_l[u]','-9.884397018','1000';'EX_met_L[u]','-2.144657123','1000';'EX_phe_L[u]','-6.083829725','1000';'EX_pro_L[u]','-11.89938505','1000';'EX_ser_L[u]','-4.424652451','1000';'EX_thr_L[u]','-3.567830759','1000';'EX_trp_L[u]','-0.685504997','1000';'EX_tyr_L[u]','-1.683306566','1000';'EX_val_L[u]','-10.37149589','1000';'EX_glc_D[u]','-27.7537245498346','1000';'EX_hco3[u]','-1.190379826','1000';'EX_phyQ[u]','-0.002172408','1000';'EX_etoh[u]','-3.237913066','1000';'EX_pheme[u]','-0.0076693','1000';'EX_oh1[u]','-0.099987502','1000';'EX_cys_L[u]','-2.846894039','1000';'EX_M02144[u]','-2.846894039','1000';'EX_h2o[u]','-55013.2623','1000';'EX_h[u]','-6.30957E-05','1000'};


%Western European Diet
dietConstraints{1}={'EX_fru[u]','-0.14899','1000';'EX_glc_D[u]','-0.14899','1000';'EX_gal[u]','-0.14899','1000';'EX_man[u]','-0.14899','1000';'EX_mnl[u]','-0.14899','1000';'EX_fuc_L[u]','-0.14899','1000';'EX_glcn[u]','-0.14899','1000';'EX_rmn[u]','-0.14899','1000';'EX_arab_L[u]','-0.17878','1000';'EX_drib[u]','-0.17878','1000';'EX_rib_D[u]','-0.17878','1000';'EX_xyl_D[u]','-0.17878','1000';'EX_oxa[u]','-0.44696','1000';'EX_lcts[u]','-0.074493','1000';'EX_malt[u]','-0.074493','1000';'EX_sucr[u]','-0.074493','1000';'EX_melib[u]','-0.074493','1000';'EX_cellb[u]','-0.074493','1000';'EX_tre[u]','-0.074493','1000';'EX_strch1[u]','-0.25734','1000';'EX_amylopect900[u]','-1.5673e-05','1000';'EX_amylose300[u]','-4.7019e-05','1000';'EX_arabinan101[u]','-0.00016628','1000';'EX_arabinogal[u]','-2.1915e-05','1000';'EX_arabinoxyl[u]','-0.00030665','1000';'EX_bglc[u]','-7.05e-08','1000';'EX_cellul[u]','-2.8212e-05','1000';'EX_dextran40[u]','-0.00017632','1000';'EX_galmannan[u]','-1.4106e-05','1000';'EX_glcmannan[u]','-3.2881e-05','1000';'EX_homogal[u]','-0.00012823','1000';'EX_inulin[u]','-0.00047019','1000';'EX_kestopt[u]','-0.0028212','1000';'EX_levan1000[u]','-1.4106e-05','1000';'EX_lmn30[u]','-0.00047019','1000';'EX_lichn[u]','-8.2976e-05','1000';'EX_pect[u]','-3.3387e-05','1000';'EX_pullulan1200[u]','-1.1755e-05','1000';'EX_raffin[u]','-0.0047019','1000';'EX_rhamnogalurI[u]','-1.4492e-05','1000';'EX_rhamnogalurII[u]','-0.00026699','1000';'EX_starch1200[u]','-1.1755e-05','1000';'EX_xylan[u]','-3.2059e-05','1000';'EX_xyluglc[u]','-1.3146e-05','1000';'EX_arachd[u]','-0.003328','1000';'EX_chsterol[u]','-0.004958','1000';'EX_glyc[u]','-1.7997','1000';'EX_hdca[u]','-0.39637','1000';'EX_hdcea[u]','-0.036517','1000';'EX_lnlc[u]','-0.35911','1000';'EX_lnlnca[u]','-0.017565','1000';'EX_lnlncg[u]','-0.017565','1000';'EX_ocdca[u]','-0.16928','1000';'EX_ocdcea[u]','-0.68144','1000';'EX_octa[u]','-0.012943','1000';'EX_ttdca[u]','-0.068676','1000';'EX_ala_L[u]','-1','1000';'EX_cys_L[u]','-1','1000';'EX_ser_L[u]','-1','1000';'EX_arg_L[u]','-0.15','1000';'EX_his_L[u]','-0.15','1000';'EX_ile_L[u]','-0.15','1000';'EX_leu_L[u]','-0.15','1000';'EX_lys_L[u]','-0.15','1000';'EX_asn_L[u]','-0.225','1000';'EX_asp_L[u]','-0.225','1000';'EX_thr_L[u]','-0.225','1000';'EX_glu_L[u]','-0.18','1000';'EX_met_L[u]','-0.18','1000';'EX_gln_L[u]','-0.18','1000';'EX_pro_L[u]','-0.18','1000';'EX_val_L[u]','-0.18','1000';'EX_phe_L[u]','-1','1000';'EX_tyr_L[u]','-1','1000';'EX_gly[u]','-0.45','1000';'EX_trp_L[u]','-0.08182','1000';'EX_12dgr180[u]','-1','1000';'EX_26dap_M[u]','-1','1000';'EX_2dmmq8[u]','-1','1000';'EX_2obut[u]','-1','1000';'EX_3mop[u]','-1','1000';'EX_4abz[u]','-1','1000';'EX_4hbz[u]','-1','1000';'EX_5aop[u]','-1','1000';'EX_ac[u]','-1','1000';'EX_acald[u]','-1','1000';'EX_acgam[u]','-1','1000';'EX_acmana[u]','-1','1000';'EX_acnam[u]','-1','1000';'EX_ade[u]','-1','1000';'EX_adn[u]','-1','1000';'EX_adocbl[u]','-1','1000';'EX_akg[u]','-1','1000';'EX_ala_D[u]','-1','1000';'EX_amet[u]','-1','1000';'EX_amp[u]','-1','1000';'EX_anth[u]','-1','1000';'EX_arab_D[u]','-1','1000';'EX_avite1[u]','-1','1000';'EX_btn[u]','-1','1000';'EX_ca2[u]','-1','1000';'EX_cbl1[u]','-1','1000';'EX_cgly[u]','-1','1000';'EX_chol[u]','-1','1000';'EX_chor[u]','-1','1000';'EX_cit[u]','-1','1000';'EX_cl[u]','-1','1000';'EX_cobalt2[u]','-1','1000';'EX_csn[u]','-1','1000';'EX_cu2[u]','-1','1000';'EX_cytd[u]','-1','1000';'EX_dad_2[u]','-1','1000';'EX_dcyt[u]','-1','1000';'EX_ddca[u]','-1','1000';'EX_dgsn[u]','-1','1000';'EX_etoh[u]','-1','1000';'EX_fald[u]','-1','1000';'EX_fe2[u]','-1','1000';'EX_fe3[u]','-1','1000';'EX_fe3dcit[u]','-1','1000';'EX_fol[u]','-1','1000';'EX_for[u]','-1','1000';'EX_fum[u]','-1','1000';'EX_gam[u]','-1','1000';'EX_glu_D[u]','-1','1000';'EX_glyc3p[u]','-1','1000';'EX_gsn[u]','-1','1000';'EX_gthox[u]','-1','1000';'EX_gthrd[u]','-1','1000';'EX_gua[u]','-1','1000';'EX_h[u]','-1','1000';'EX_h2[u]','-1','1000';'EX_h2s[u]','-1','1000';'EX_hom_L[u]','-1','1000';'EX_hxan[u]','-1','1000';'EX_indole[u]','-1','1000';'EX_ins[u]','-1','1000';'EX_k[u]','-1','1000';'EX_lac_L[u]','-1','1000';'EX_lanost[u]','-1','1000';'EX_mal_L[u]','-1','1000';'EX_metsox_S_L[u]','-1','1000';'EX_mg2[u]','-1','1000';'EX_mn2[u]','-1','1000';'EX_mobd[u]','-1','1000';'EX_mqn7[u]','-1','1000';'EX_mqn8[u]','-1','1000';'EX_na1[u]','-1','1000';'EX_nac[u]','-1','1000';'EX_ncam[u]','-1','1000';'EX_nmn[u]','-1','1000';'EX_no2[u]','-1','1000';'EX_no2[u]','-1','1000';'EX_no3[u]','-1','1000';'EX_orn[u]','-1','1000';'EX_pheme[u]','-1','1000';'EX_phyQ[u]','-1','1000';'EX_pi[u]','-1','1000';'EX_pime[u]','-1','1000';'EX_pnto_R[u]','-1','1000';'EX_ptrc[u]','-1','1000';'EX_pydam[u]','-1','1000';'EX_pydx[u]','-1','1000';'EX_pydx5p[u]','-1','1000';'EX_pydxn[u]','-1','1000';'EX_q8[u]','-1','1000';'EX_retinol[u]','-1','1000';'EX_ribflv[u]','-1','1000';'EX_sel[u]','-1','1000';'EX_sheme[u]','-1','1000';'EX_so4[u]','-1','1000';'EX_spmd[u]','-1','1000';'EX_succ[u]','-1','1000';'EX_thf[u]','-1','1000';'EX_thm[u]','-1','1000';'EX_thymd[u]','-1','1000';'EX_ura[u]','-1','1000';'EX_uri[u]','-1','1000';'EX_vitd3[u]','-1','1000';'EX_xan[u]','-1','1000';'EX_zn2[u]','-1','1000';'EX_meoh[u]','-10','1000';'EX_h2o[u]','-10','1000'};

%YCAG  Human2
%dietConstraints{3}={'EX_pi[u]','-7.605312000000001','1000';'EX_k[u]','-15.22291','1000';'EX_4abz[u]','-0.00021899999999999998','1000';'EX_arab_L[u]','-15.144173','1000';'EX_arg_L[u]','-1.9804739999999998','1000';'EX_arg_D[u]','-1.9804739999999998','1000';'EX_btn[u]','-5.47e-05','1000';'EX_ca[u]','-0.8109569999999999','1000';'EX_cl[u]','-24.457700999999997','1000';'EX_cbl1[u]','-7.38e-06','1000';'EX_cys_L[u]','-6.716115000000001','1000';'EX_h[u]','-6.344423','1000';'EX_so4[u]','-3.8777579999999996','1000';'EX_h2o[u]','-27754.6489','1000';'EX_fol[u]','-0.00011300000000000001','1000';'EX_na1[u]','-78.695299','1000';'EX_glc_D[u]','-0.0','1000';'EX_gal[u]','-15.144173','1000';'EX_mg2[u]','-2.243342','1000';'EX_mn2[u]','-0.01','1000';'EX_hco3[u]','-47.613378999999995','1000';'EX_nh4[u]','-1.8318919999999999','1000';'EX_fe3[u]','-0.0218','1000';'EX_pydam[u]','-0.000622','1000';'EX_thm[u]','-0.00104','1000';'EX_ala_D[u]','-4.377525','1000';'EX_ala_L[u]','-4.377525','1000';'EX_asp_L[u]','-5.327528','1000';'EX_gln_L[u]','-0.17105499999999998','1000';'EX_gln_D[u]','-0.17105499999999998','1000';'EX_glu_D[u]','-6.353279','1000';'EX_glu_L[u]','-6.353279','1000';'EX_gly[u]','-3.796496','1000';'EX_his_L[u]','-0.918456','1000';'EX_ile_L[u]','-2.706389','1000';'EX_leu_L[u]','-4.059588000000001','1000';'EX_lys_L[u]','-2.770368','1000';'EX_met_L[u]','-0.469139','1000';'EX_phe_L[u]','-2.270109','1000';'EX_pro_L[u]','-2.1714510000000002','1000';'EX_ser_L[u]','-1.807963','1000';'EX_asn_L[u]','-1.112662','1000';'EX_trp_L[u]','-0.1591','1000';'EX_tyr_L[u]','-0.88305','1000';'EX_val_L[u]','-3.0516560000000004','1000';'EX_ca2[u]','-0.145336','1000';'EX_fe2[u]','-0.0218','1000';'EX_thr_L[u]','-0.33579600000000004','1000';'EX_cd2[u]','-3.34e-05','1000';'EX_cobalt2[u]','-0.01','1000';'EX_cu2[u]','-0.01','1000';'EX_ni2[u]','-0.000272','1000';'EX_zn2[u]','-0.141403','1000';'EX_ade[u]','-0.163546','1000';'EX_gua[u]','-0.15301499999999998','1000';'EX_csn[u]','-0.09809999999999999','1000';'EX_ura[u]','-0.140738','1000';'EX_man[u]','-0.63556','1000';'EX_pnto_R[u]','-0.0012100000000000001','1000';'EX_ascb_L[u]','-4.26e-05','1000';'EX_nac[u]','-0.00816','1000';'EX_pydx[u]','-0.000222','1000';'EX_chol[u]','-0.007679999999999999','1000';'EX_adocbl[u]','-1.2899999999999999e-09','1000';'EX_mobd[u]','-0.000154','1000';'EX_n2[u]','-800.0','1000';'EX_h2[u]','-100.0','1000';'EX_26dap_M[u]','-10e-4','1000';'EX_h2s[u]','-10e-4','1000';'EX_ribflv[u]','-10e-4','1000';'EX_thymd[u]','-10e-4','1000';'EX_12dgr180[u]','-10e-4','1000';'EX_hxan[u]','-10e-4','1000';'EX_ddca[u]','-10e-4','1000';'EX_ttdca[u]','-10e-4','1000';'EX_pheme[u]','-10e-4','1000';'EX_q8[u]','-10e-4','1000';'EX_arab_D[u]','-10e-4','1000';'EX_cgly[u]','-10e-4','1000';'EX_cit[u]','-10e-4','1000';'EX_no2[u]','-10e-4','1000';'EX_no3[u]','-10e-4','1000';'EX_ocdca[u]','-10e-4','1000';'EX_orn[u]','-10e-4','1000';'EX_ptrc[u]','-10e-4','1000';'EX_sheme[u]','-10e-4','1000';'EX_spmd[u]','-10e-4','1000';'EX_mqn7[u]','-10e-4','1000'};

%YCGMS Human3
%dietConstraints{3}={'EX_pi[u]','-7.605312000000001','1000';'EX_k[u]','-15.22291','1000';'EX_4abz[u]','-0.00021899999999999998','1000';'EX_arg_L[u]','-1.9804739999999998','1000';'EX_arg_D[u]','-1.9804739999999998','1000';'EX_btn[u]','-5.47e-05','1000';'EX_ca[u]','-0.8109569999999999','1000';'EX_cl[u]','-24.457700999999997','1000';'EX_cbl1[u]','-7.38e-06','1000';'EX_cellb[u]','-0.00584','1000';'EX_cys_L[u]','-6.716115000000001','1000';'EX_h[u]','-6.344423','1000';'EX_so4[u]','-3.8777579999999996','1000';'EX_h2o[u]','-27754.6489','1000';'EX_fol[u]','-0.00011300000000000001','1000';'EX_na1[u]','-78.695299','1000';'EX_glc_D[u]','-0.0111','1000';'EX_malt[u]','-0.00555','1000';'EX_mg2[u]','-2.243342','1000';'EX_mn2[u]','-0.01','1000';'EX_hco3[u]','-47.613378999999995','1000';'EX_nh4[u]','-1.8318919999999999','1000';'EX_fe3[u]','-0.0218','1000';'EX_pydam[u]','-0.000622','1000';'EX_ribflv[u]','-0.000133','1000';'EX_strch2[u]','-0.0010092','1000';'EX_thm[u]','-0.0011489999999999998','1000';'EX_ala_D[u]','-4.377525','1000';'EX_ala_L[u]','-4.377525','1000';'EX_asp_L[u]','-5.327528','1000';'EX_gln_L[u]','-0.17105499999999998','1000';'EX_gln_D[u]','-0.17105499999999998','1000';'EX_glu_D[u]','-6.353279','1000';'EX_glu_L[u]','-6.353279','1000';'EX_gly[u]','-3.796496','1000';'EX_his_L[u]','-0.918456','1000';'EX_ile_L[u]','-2.706389','1000';'EX_leu_L[u]','-4.059588000000001','1000';'EX_lys_L[u]','-2.770368','1000';'EX_met_L[u]','-0.469139','1000';'EX_phe_L[u]','-2.270109','1000';'EX_pro_L[u]','-2.1714510000000002','1000';'EX_ser_L[u]','-1.807963','1000';'EX_asn_L[u]','-1.112662','1000';'EX_trp_L[u]','-0.1591','1000';'EX_tyr_L[u]','-0.88305','1000';'EX_val_L[u]','-3.0516560000000004','1000';'EX_ca2[u]','-0.145336','1000';'EX_fe2[u]','-0.0218','1000';'EX_thr_L[u]','-0.33579600000000004','1000';'EX_cd2[u]','-3.34e-05','1000';'EX_cobalt2[u]','-0.01','1000';'EX_cu2[u]','-0.01','1000';'EX_ni2[u]','-0.000272','1000';'EX_zn2[u]','-0.141403','1000';'EX_ade[u]','-0.163546','1000';'EX_gua[u]','-0.15301499999999998','1000';'EX_csn[u]','-0.09809999999999999','1000';'EX_ura[u]','-0.140738','1000';'EX_man[u]','-0.63556','1000';'EX_pnto_R[u]','-0.0012100000000000001','1000';'EX_ascb_L[u]','-4.26e-05','1000';'EX_nac[u]','-0.00816','1000';'EX_pydx[u]','-0.000222','1000';'EX_chol[u]','-0.007679999999999999','1000';'EX_adocbl[u]','-1.2899999999999999e-09','1000';'EX_mobd[u]','-0.000154','1000';'EX_n2[u]','-800.0','1000';'EX_h2[u]','-100.0','1000';'EX_26dap_M[u]','-10e-4','1000';'EX_h2s[u]','-10e-4','1000';'EX_thymd[u]','-10e-4','1000';'EX_12dgr180[u]','-10e-4','1000';'EX_hxan[u]','-10e-4','1000';'EX_ddca[u]','-10e-4','1000';'EX_ttdca[u]','-10e-4','1000';'EX_pheme[u]','-10e-4','1000';'EX_q8[u]','-10e-4','1000';'EX_arab_D[u]','-10e-4','1000';'EX_cgly[u]','-10e-4','1000';'EX_cit[u]','-10e-4','1000';'EX_no2[u]','-10e-4','1000';'EX_no3[u]','-10e-4','1000';'EX_ocdca[u]','-10e-4','1000';'EX_orn[u]','-10e-4','1000';'EX_ptrc[u]','-10e-4','1000';'EX_sheme[u]','-10e-4','1000';'EX_spmd[u]','-10e-4','1000';'EX_mqn7[u]','-10e-4','1000'};

%YCFA  Human4
%dietConstraints{4}={'EX_pi[u]','-7.605312000000001','1000';'EX_k[u]','-15.22291','1000';'EX_4abz[u]','-0.00021899999999999998','1000';'EX_ac[u]','-33.0','1000';'EX_arg_L[u]','-1.9804739999999998','1000';'EX_arg_D[u]','-1.9804739999999998','1000';'EX_btn[u]','-5.47e-05','1000';'EX_ca[u]','-0.8109569999999999','1000';'EX_cl[u]','-24.457700999999997','1000';'EX_cbl1[u]','-7.38e-06','1000';'EX_cys_L[u]','-6.716115000000001','1000';'EX_h[u]','-6.344423','1000';'EX_so4[u]','-3.8777579999999996','1000';'EX_h2o[u]','-27754.6489','1000';'EX_fol[u]','-0.00011300000000000001','1000';'EX_na1[u]','-78.695299','1000';'EX_glc_D[u]','-27.753107999999997','1000';'EX_mg2[u]','-2.243342','1000';'EX_mn2[u]','-0.01','1000';'EX_hco3[u]','-47.613378999999995','1000';'EX_nh4[u]','-1.8318919999999999','1000';'EX_fe3[u]','-0.0218','1000';'EX_pydam[u]','-0.000622','1000';'EX_thm[u]','-0.00104','1000';'EX_ppa[u]','-9.0','1000';'EX_isobut[u]','-1.0','1000';'EX_isoval[u]','-1.0','1000';'EX_M03134[u]','-1.0','1000';'EX_ala_D[u]','-4.377525','1000';'EX_ala_L[u]','-4.377525','1000';'EX_asp_L[u]','-5.327528','1000';'EX_gln_L[u]','-0.17105499999999998','1000';'EX_gln_D[u]','-0.17105499999999998','1000';'EX_glu_D[u]','-6.353279','1000';'EX_glu_L[u]','-6.353279','1000';'EX_gly[u]','-3.796496','1000';'EX_his_L[u]','-0.918456','1000';'EX_ile_L[u]','-2.706389','1000';'EX_leu_L[u]','-4.059588000000001','1000';'EX_lys_L[u]','-2.770368','1000';'EX_met_L[u]','-0.469139','1000';'EX_phe_L[u]','-2.270109','1000';'EX_pro_L[u]','-2.1714510000000002','1000';'EX_ser_L[u]','-1.807963','1000';'EX_asn_L[u]','-1.112662','1000';'EX_trp_L[u]','-0.1591','1000';'EX_tyr_L[u]','-0.88305','1000';'EX_val_L[u]','-3.0516560000000004','1000';'EX_ca2[u]','-0.145336','1000';'EX_fe2[u]','-0.0218','1000';'EX_thr_L[u]','-0.33579600000000004','1000';'EX_cd2[u]','-3.34e-05','1000';'EX_cobalt2[u]','-0.01','1000';'EX_cu2[u]','-0.01','1000';'EX_ni2[u]','-0.000272','1000';'EX_zn2[u]','-0.141403','1000';'EX_ade[u]','-0.163546','1000';'EX_gua[u]','-0.15301499999999998','1000';'EX_csn[u]','-0.09809999999999999','1000';'EX_ura[u]','-0.140738','1000';'EX_man[u]','-0.63556','1000';'EX_pnto_R[u]','-0.0012100000000000001','1000';'EX_ascb_L[u]','-4.26e-05','1000';'EX_nac[u]','-0.00816','1000';'EX_pydx[u]','-0.000222','1000';'EX_chol[u]','-0.007679999999999999','1000';'EX_adocbl[u]','-1.2899999999999999e-09','1000';'EX_mobd[u]','-0.000154','1000';'EX_n2[u]','-800.0','1000';'EX_h2[u]','-100.0','1000';'EX_26dap_M[u]','-10e-4','1000';'EX_h2s[u]','-10e-4','1000';'EX_ribflv[u]','-10e-4','1000';'EX_thymd[u]','-10e-4','1000';'EX_12dgr180[u]','-10e-4','1000';'EX_hxan[u]','-10e-4','1000';'EX_ddca[u]','-10e-4','1000';'EX_ttdca[u]','-10e-4','1000';'EX_pheme[u]','-10e-4','1000';'EX_q8[u]','-10e-4','1000';'EX_arab_D[u]','-10e-4','1000';'EX_cgly[u]','-10e-4','1000';'EX_cit[u]','-10e-4','1000';'EX_no2[u]','-10e-4','1000';'EX_no3[u]','-10e-4','1000';'EX_ocdca[u]','-10e-4','1000';'EX_orn[u]','-10e-4','1000';'EX_ptrc[u]','-10e-4','1000';'EX_sheme[u]','-10e-4','1000';'EX_spmd[u]','-10e-4','1000';'EX_mqn7[u]','-10e-4','1000'};

%mMCB Human5 Human6
%dietConstraints{5}={'EX_pi[u]','-1507.9056229999999','1000';'EX_k[u]','-16.288865','1000';'EX_ac[u]','-500.0','1000';'EX_arg_L[u]','-2.296201','1000';'EX_arg_D[u]','-2.296201','1000';'EX_btn[u]','-1.6499999999999998e-05','1000';'EX_ca[u]','-0.8109569999999999','1000';'EX_cl[u]','-82.428343','1000';'EX_cbl1[u]','-0.0','1000';'EX_co2[u]','-100.0','1000';'EX_cys_L[u]','-2.967312','1000';'EX_h[u]','-4502.5378089999995','1000';'EX_fe[u]','-0.018000000000000002','1000';'EX_so4[u]','-3.970123','1000';'EX_h2o[u]','-33305.899594','1000';'EX_for[u]','-350.0','1000';'EX_na1[u]','-1947.443099','1000';'EX_fru[u]','-50.0','1000';'EX_glc_D[u]','-0.0','1000';'EX_h2[u]','-100.0','1000';'EX_mg2[u]','-2.509357','1000';'EX_mn2[u]','-0.334006','1000';'EX_mqn8[u]','-0.029','1000';'EX_n2[u]','-800.0','1000';'EX_hco3[u]','-2.380669','1000';'EX_oh1[u]','-1500.0','1000';'EX_nh4[u]','-1.380955','1000';'EX_fe3[u]','-0.0261','1000';'EX_pydam[u]','-0.0','1000';'EX_thm[u]','-0.00124','1000';'EX_ala_D[u]','-5.112723','1000';'EX_ala_L[u]','-5.112723','1000';'EX_asp_L[u]','-6.176431','1000';'EX_gln_L[u]','-0.198483','1000';'EX_gln_D[u]','-0.198483','1000';'EX_glu_D[u]','-7.386132','1000';'EX_glu_L[u]','-7.386132','1000';'EX_gly[u]','-4.415926000000001','1000';'EX_his_L[u]','-1.0666980000000001','1000';'EX_ile_L[u]','-3.140936','1000';'EX_leu_L[u]','-4.707596','1000';'EX_lys_L[u]','-3.2252549999999998','1000';'EX_met_L[u]','-0.546212','1000';'EX_phe_L[u]','-2.630299','1000';'EX_pro_L[u]','-2.5188829999999998','1000';'EX_ser_L[u]','-2.098188','1000';'EX_asn_L[u]','-1.289022','1000';'EX_trp_L[u]','-0.186018','1000';'EX_tyr_L[u]','-1.023786','1000';'EX_val_L[u]','-3.54675','1000';'EX_ca2[u]','-0.167539','1000';'EX_fe2[u]','-0.0261','1000';'EX_thr_L[u]','-0.40295499999999995','1000';'EX_cd2[u]','-4e-05','1000';'EX_cobalt2[u]','-0.01','1000';'EX_cu2[u]','-0.01','1000';'EX_ni2[u]','-0.000327','1000';'EX_zn2[u]','-0.19584','1000';'EX_ade[u]','-0.19625499999999999','1000';'EX_gua[u]','-0.183618','1000';'EX_csn[u]','-0.117728','1000';'EX_ura[u]','-0.16888499999999998','1000';'EX_man[u]','-0.762672','1000';'EX_pnto_R[u]','-0.0014500000000000001','1000';'EX_ascb_L[u]','-5.11e-05','1000';'EX_nac[u]','-0.0098','1000';'EX_pydx[u]','-0.000266','1000';'EX_chol[u]','-0.009219999999999999','1000';'EX_adocbl[u]','-1.55e-09','1000';'EX_mobd[u]','-0.000184','1000';'EX_h2[u]','-100.0','1000';'EX_26dap_M[u]','-10e-4','1000';'EX_h2s[u]','-10e-4','1000';'EX_ribflv[u]','-10e-4','1000';'EX_thymd[u]','-10e-4','1000';'EX_fol[u]','-10e-4','1000';'EX_12dgr180[u]','-10e-4','1000';'EX_hxan[u]','-10e-4','1000';'EX_ddca[u]','-10e-4','1000';'EX_ttdca[u]','-10e-4','1000';'EX_pheme[u]','-10e-4','1000';'EX_q8[u]','-10e-4','1000';'EX_arab_D[u]','-10e-4','1000';'EX_cgly[u]','-10e-4','1000';'EX_cit[u]','-10e-4','1000';'EX_no2[u]','-10e-4','1000';'EX_no3[u]','-10e-4','1000';'EX_ocdca[u]','-10e-4','1000';'EX_orn[u]','-10e-4','1000';'EX_ptrc[u]','-10e-4','1000';'EX_sheme[u]','-10e-4','1000';'EX_spmd[u]','-0e-4','1000';'EX_mqn7[u]','-10e-4','1000';'EX_cbl1[u]','-10e-4','1000'};

%%
% need to determine reduced cost for PYG media 





[solution_opt] = optimizeCbModel(model);





 %% 
% NOTE: if you design your own diet, make sure that exchange reaction abbreviations 
% correspond to the lumen exchanges in the joint models ('EX_compound[u]').
% 
% Define what counts as significant difference between single growth of the 
% microbes and growth when joined with another microbe-here we choose 10%.

sigD = 0.1;
%% 
% Simulate the pairwise interactions on the four dietary conditions.

for i = 1:length(conditions)
    % assign dietary constraints
    [pairwiseInteractions,pairwiseSolutions]=simulatePairwiseInteractions(resPath,'inputDiet',dietConstraints{i},'sigD',sigD,'saveSolutionsFlag', true,'numWorkers', numWorkers);
Interactions.(conditions{i})=pairwiseInteractions;
end
%% Analysis of computed pairwise interactions
% The computed microbe-microbe interactions will be plotted by type and analyzed 
% in the context of the taxonomy of the joined strains. There are six possible 
% types of interactions total that can result in increased growth (+), no change 
% in growth (=) or decreased growth (-) compared with the single condition for 
% each joined microbe.
%% 
% * Competition (-/-)
% * Parasitism (+/-)
% * Amensalism (=/-)
% * Neutralism (=/=)
% * Commensalism (+/=)
% * Mutualism (+/+)
%% 
% This results in nine different outcomes total from the perspective of each 
% joined microbe.
% 
% Plot the percentage of interactions computed.

figure('rend','painters','pos',[10 10 900 600])
typesIA=unique(pairwiseInteractions(2:end,10));
for i = 1:length(conditions)
    pairwiseInteractions=Interactions.(conditions{i});
    listIA=pairwiseInteractions(2:end,10);
    for j=1:length(typesIA)
        dat(j)=sum(strcmp(listIA(:),typesIA{j}));
    end
    subplot(2,2,i)
    pie(dat)
    set(gca,'FontSize',10)
    h=title(conditions{i});
    set(h,'interpreter','none')
    title(conditions{i})
end
legend1=legend(typesIA);
set(legend1,'Position',[0.42 0.45 0.2 0.2],'FontSize',12)
sgtitle('Percentage of computed pairwise interactions')
%% 
% Next, the percentage of interactions will be calculated on different taxon 
% levels (genus, family, order, class, phylum) using the taxon information contained 
% in AGORA_infoFile.xlsx. Here, the interactions will be considered from the perspective 
% of each joined microbe resulting in nine possible interactions total.
% 
% Calculate the percentage of interactions predicted for each taxon included 
% in the list of microbes analyzed.

for i = 1:length(conditions)
    pairwiseInteractions=Interactions.(conditions{i});
    [InteractionsByTaxon]=calculateInteractionsByTaxon(pairwiseInteractions,infoFile);
    TaxonSummaries.(conditions{i})=InteractionsByTaxon;
end
%% 
% Combine the four conditions into one structure.

InteractionsByTaxonCombined=struct;
for i = 1:length(conditions)
    InteractionsByTaxon=TaxonSummaries.(conditions{i});
    taxLevels=fieldnames(InteractionsByTaxon);
    if i==1
        for j=1:length(taxLevels)
            InteractionsByTaxonCombined.(taxLevels{j})=InteractionsByTaxon.(taxLevels{j});
            InteractionsByTaxonCombined.(taxLevels{j})(2:end,1)=strcat(InteractionsByTaxonCombined.(taxLevels{j})(2:end,1),'_',conditions{i});
        end
    else
        for j=1:length(taxLevels)
            rowLength=size(InteractionsByTaxonCombined.(taxLevels{j}),1);
            InteractionsByTaxonCombined.(taxLevels{j})=[InteractionsByTaxonCombined.(taxLevels{j});InteractionsByTaxon.(taxLevels{j})(2:end,:)];
            InteractionsByTaxonCombined.(taxLevels{j})(rowLength+1:end,1)=strcat(InteractionsByTaxonCombined.(taxLevels{j})(rowLength+1:end,1),'_',conditions{i});
        end
    end
end
%% 
% Let us plot the distributions of interactions for all dietary conditions combined 
% on the level of genera as an example. Note: The xticklabels/yticklabels function 
% is only available in MATLAB R2016b or newer. Older versions of MATLAB will be 
% unable to display the labels.

for i=5
    xlabels=InteractionsByTaxonCombined.(taxLevels{i})(1,2:end);
    ylabels=InteractionsByTaxonCombined.(taxLevels{i})(2:end,1);
    data=string(InteractionsByTaxonCombined.(taxLevels{i})(2:end,2:end));
    data=str2double(data);
    figure;
    imagesc(data)
    colormap('hot')
    colorbar
    set(gca,'xtick',1:length(xlabels));
    xticklabels(xlabels);
    set(gca,'ytick',1:length(ylabels));
    yticklabels(ylabels);
    xtickangle(90)
    set(gca,'TickLabelInterpreter', 'none');
    title(taxLevels{i})
end
%% Pareto optimality analysis
% Another way to analyze the metabolic interactions between two microbes in 
% Pareto optimality analysis. In this method, the tradeoff between two competing 
% objectives (e.g., the biomasses of two joined microbes) is calculated. The resulting 
% Pareto frontier depicts all possible outcomes of co-growth between the two microbes 
% under the given constraints.
% 
% Let us compute the Pareto frontier for five randomly chosen pairs from the 
% list of AGORA models.

modelInd = randi([2 length(infoFile)],2,5);
%% 
% The Pareto frontier will be computed on the Western diet without oxygen.

%dietConstraints{1}={'EX_fru[u]','-0.14899','1000';'EX_glc_D[u]','-0.14899','1000';'EX_gal[u]','-0.14899','1000';'EX_man[u]','-0.14899','1000';'EX_mnl[u]','-0.14899','1000';'EX_fuc_L[u]','-0.14899','1000';'EX_glcn[u]','-0.14899','1000';'EX_rmn[u]','-0.14899','1000';'EX_arab_L[u]','-0.17878','1000';'EX_drib[u]','-0.17878','1000';'EX_rib_D[u]','-0.17878','1000';'EX_xyl_D[u]','-0.17878','1000';'EX_oxa[u]','-0.44696','1000';'EX_lcts[u]','-0.074493','1000';'EX_malt[u]','-0.074493','1000';'EX_sucr[u]','-0.074493','1000';'EX_melib[u]','-0.074493','1000';'EX_cellb[u]','-0.074493','1000';'EX_tre[u]','-0.074493','1000';'EX_strch1[u]','-0.25734','1000';'EX_amylopect900[u]','-1.5673e-05','1000';'EX_amylose300[u]','-4.7019e-05','1000';'EX_arabinan101[u]','-0.00016628','1000';'EX_arabinogal[u]','-2.1915e-05','1000';'EX_arabinoxyl[u]','-0.00030665','1000';'EX_bglc[u]','-7.05e-08','1000';'EX_cellul[u]','-2.8212e-05','1000';'EX_dextran40[u]','-0.00017632','1000';'EX_galmannan[u]','-1.4106e-05','1000';'EX_glcmannan[u]','-3.2881e-05','1000';'EX_homogal[u]','-0.00012823','1000';'EX_inulin[u]','-0.00047019','1000';'EX_kestopt[u]','-0.0028212','1000';'EX_levan1000[u]','-1.4106e-05','1000';'EX_lmn30[u]','-0.00047019','1000';'EX_lichn[u]','-8.2976e-05','1000';'EX_pect[u]','-3.3387e-05','1000';'EX_pullulan1200[u]','-1.1755e-05','1000';'EX_raffin[u]','-0.0047019','1000';'EX_rhamnogalurI[u]','-1.4492e-05','1000';'EX_rhamnogalurII[u]','-0.00026699','1000';'EX_starch1200[u]','-1.1755e-05','1000';'EX_xylan[u]','-3.2059e-05','1000';'EX_xyluglc[u]','-1.3146e-05','1000';'EX_arachd[u]','-0.003328','1000';'EX_chsterol[u]','-0.004958','1000';'EX_glyc[u]','-1.7997','1000';'EX_hdca[u]','-0.39637','1000';'EX_hdcea[u]','-0.036517','1000';'EX_lnlc[u]','-0.35911','1000';'EX_lnlnca[u]','-0.017565','1000';'EX_lnlncg[u]','-0.017565','1000';'EX_ocdca[u]','-0.16928','1000';'EX_ocdcea[u]','-0.68144','1000';'EX_octa[u]','-0.012943','1000';'EX_ttdca[u]','-0.068676','1000';'EX_ala_L[u]','-1','1000';'EX_cys_L[u]','-1','1000';'EX_ser_L[u]','-1','1000';'EX_arg_L[u]','-0.15','1000';'EX_his_L[u]','-0.15','1000';'EX_ile_L[u]','-0.15','1000';'EX_leu_L[u]','-0.15','1000';'EX_lys_L[u]','-0.15','1000';'EX_asn_L[u]','-0.225','1000';'EX_asp_L[u]','-0.225','1000';'EX_thr_L[u]','-0.225','1000';'EX_glu_L[u]','-0.18','1000';'EX_met_L[u]','-0.18','1000';'EX_gln_L[u]','-0.18','1000';'EX_pro_L[u]','-0.18','1000';'EX_val_L[u]','-0.18','1000';'EX_phe_L[u]','-1','1000';'EX_tyr_L[u]','-1','1000';'EX_gly[u]','-0.45','1000';'EX_trp_L[u]','-0.08182','1000';'EX_12dgr180[u]','-1','1000';'EX_26dap_M[u]','-1','1000';'EX_2dmmq8[u]','-1','1000';'EX_2obut[u]','-1','1000';'EX_3mop[u]','-1','1000';'EX_4abz[u]','-1','1000';'EX_4hbz[u]','-1','1000';'EX_5aop[u]','-1','1000';'EX_ac[u]','-1','1000';'EX_acald[u]','-1','1000';'EX_acgam[u]','-1','1000';'EX_acmana[u]','-1','1000';'EX_acnam[u]','-1','1000';'EX_ade[u]','-1','1000';'EX_adn[u]','-1','1000';'EX_adocbl[u]','-1','1000';'EX_akg[u]','-1','1000';'EX_ala_D[u]','-1','1000';'EX_amet[u]','-1','1000';'EX_amp[u]','-1','1000';'EX_anth[u]','-1','1000';'EX_arab_D[u]','-1','1000';'EX_avite1[u]','-1','1000';'EX_btn[u]','-1','1000';'EX_ca2[u]','-1','1000';'EX_cbl1[u]','-1','1000';'EX_cgly[u]','-1','1000';'EX_chol[u]','-1','1000';'EX_chor[u]','-1','1000';'EX_cit[u]','-1','1000';'EX_cl[u]','-1','1000';'EX_cobalt2[u]','-1','1000';'EX_csn[u]','-1','1000';'EX_cu2[u]','-1','1000';'EX_cytd[u]','-1','1000';'EX_dad_2[u]','-1','1000';'EX_dcyt[u]','-1','1000';'EX_ddca[u]','-1','1000';'EX_dgsn[u]','-1','1000';'EX_etoh[u]','-1','1000';'EX_fald[u]','-1','1000';'EX_fe2[u]','-1','1000';'EX_fe3[u]','-1','1000';'EX_fe3dcit[u]','-1','1000';'EX_fol[u]','-1','1000';'EX_for[u]','-1','1000';'EX_fum[u]','-1','1000';'EX_gam[u]','-1','1000';'EX_glu_D[u]','-1','1000';'EX_glyc3p[u]','-1','1000';'EX_gsn[u]','-1','1000';'EX_gthox[u]','-1','1000';'EX_gthrd[u]','-1','1000';'EX_gua[u]','-1','1000';'EX_h[u]','-1','1000';'EX_h2[u]','-1','1000';'EX_h2s[u]','-1','1000';'EX_hom_L[u]','-1','1000';'EX_hxan[u]','-1','1000';'EX_indole[u]','-1','1000';'EX_ins[u]','-1','1000';'EX_k[u]','-1','1000';'EX_lac_L[u]','-1','1000';'EX_lanost[u]','-1','1000';'EX_mal_L[u]','-1','1000';'EX_metsox_S_L[u]','-1','1000';'EX_mg2[u]','-1','1000';'EX_mn2[u]','-1','1000';'EX_mobd[u]','-1','1000';'EX_mqn7[u]','-1','1000';'EX_mqn8[u]','-1','1000';'EX_na1[u]','-1','1000';'EX_nac[u]','-1','1000';'EX_ncam[u]','-1','1000';'EX_nmn[u]','-1','1000';'EX_no2[u]','-1','1000';'EX_no2[u]','-1','1000';'EX_no3[u]','-1','1000';'EX_orn[u]','-1','1000';'EX_pheme[u]','-1','1000';'EX_phyQ[u]','-1','1000';'EX_pi[u]','-1','1000';'EX_pime[u]','-1','1000';'EX_pnto_R[u]','-1','1000';'EX_ptrc[u]','-1','1000';'EX_pydam[u]','-1','1000';'EX_pydx[u]','-1','1000';'EX_pydx5p[u]','-1','1000';'EX_pydxn[u]','-1','1000';'EX_q8[u]','-1','1000';'EX_retinol[u]','-1','1000';'EX_ribflv[u]','-1','1000';'EX_sel[u]','-1','1000';'EX_sheme[u]','-1','1000';'EX_so4[u]','-1','1000';'EX_spmd[u]','-1','1000';'EX_succ[u]','-1','1000';'EX_thf[u]','-1','1000';'EX_thm[u]','-1','1000';'EX_thymd[u]','-1','1000';'EX_ura[u]','-1','1000';'EX_uri[u]','-1','1000';'EX_vitd3[u]','-1','1000';'EX_xan[u]','-1','1000';'EX_zn2[u]','-1','1000';'EX_meoh[u]','-10','1000';'EX_h2o[u]','-10','1000'};

dietConstraints{1} = {'EX_q8[u]','-1','1000';'EX_pydx[u]','-1','1000';'EX_nac[u]','-1','1000';'EX_lys_L[u]','-1','1000';'EX_hxan[u]','-1','1000';'EX_ade[u]','-1','1000';'EX_thymd[u]','-1','1000';'EX_thm[u]','-1','1000';'EX_ribflv[u]','-1','1000';'EX_pnto_R[u]','-1','1000';'EX_nac[u]','-1','1000';'EX_fol[u]','-1','1000';'EX_zn2[u]','-1','1000';'EX_mn2[u]','-1','1000';'EX_fe3[u]','-1','1000';'EX_cu2[u]','-1','1000';'EX_cobalt2[u]','-1','1000';'EX_n2[u]','-227.8571429','1000';'EX_na1[u]','-20.51163798','1000';'EX_cl[u]','-5.941065976','1000';'EX_ca2[u]','-0.173941043','1000';'EX_fe2[u]','-0.016053362','1000';'EX_mg2[u]','-0.474477191','1000';'EX_k[u]','-35.39748582','1000';'EX_so4[u]','-0.710983412','1000';'EX_pi[u]','-18.29826648','1000';'EX_ala_L[u]','-16.89227108','1000';'EX_arg_L[u]','-5.338568575','1000';'EX_asn_L[u]','-1.286718791','1000';'EX_asp_L[u]','-10.81893313','1000';'EX_Lcystin[u]','-0.187272152','1000';'EX_glu_L[u]','-18.56754922','1000';'EX_gln_L[u]','-0.205274178','1000';'EX_gly[u]','-20.3151851','1000';'EX_his_L[u]','-3.319218598','1000';'EX_ile_L[u]','-8.195159139','1000';'EX_leu_L[u]','-11.2826377','1000';'EX_lys_l[u]','-9.884397018','1000';'EX_met_L[u]','-2.144657123','1000';'EX_phe_L[u]','-6.083829725','1000';'EX_pro_L[u]','-11.89938505','1000';'EX_ser_L[u]','-4.424652451','1000';'EX_thr_L[u]','-3.567830759','1000';'EX_trp_L[u]','-0.685504997','1000';'EX_tyr_L[u]','-1.683306566','1000';'EX_val_L[u]','-10.37149589','1000';'EX_glc_D[u]','-27.7537245498346','1000';'EX_hco3[u]','-1.190379826','1000';'EX_phyQ[u]','-0.002172408','1000';'EX_etoh[u]','-3.237913066','1000';'EX_pheme[u]','-0.0076693','1000';'EX_oh1[u]','-0.099987502','1000';'EX_cys_L[u]','-2.846894039','1000';'EX_M02144[u]','-2.846894039','1000';'EX_h2o[u]','-55013.2623','1000';'EX_h[u]','-6.30957E-05','1000'};



%% 
% By default, the points of the frontier will be generated at steps of 0.001.

dinc=0.001;
%% 
% Perform the Pareto optimality analysis for the five pairs. The shape of the 
% computed Pareto frontier, which represents all possible optimal solutions of 
% simultaneously optimized growth, depends on the metabolic networks of the two 
% joined microbes.

for i=1:size(modelInd,2)
    models={};
    model=readCbModel([modPath filesep infoFile{modelInd(1,i),1} '.mat']);
    models{1,1}=model;
    bioID{1,1}=model.rxns{find(strncmp(model.rxns, 'bio', 3)),1}; % made modification in this line to grab biomass from agora2.mat files, because not always defined as just biomass
    nameTagsModels{1,1}=strcat(infoFile{modelInd(1,i),1},'_');
    model=readCbModel([modPath filesep infoFile{modelInd(2,i),1} '.mat']);
    models{2,1}=model;
    nameTagsModels{2,1}=strcat(infoFile{modelInd(2,i),1},'_');
    bioID{2,1}=model.rxns{find(strncmp(model.rxns, 'bio', 3)),1}; % made modification in this line to grab biomass from agora2.mat files, because not always defined as just biomass
    [pairedModel] = createMultipleSpeciesModel(models,'nameTagsModels',nameTagsModels);
    [pairedModel]=coupleRxnList2Rxn(pairedModel,pairedModel.rxns(strmatch(nameTagsModels{1,1},pairedModel.rxns)),strcat(infoFile{modelInd(1,i),1},'_',bioID{1,1}));
    [pairedModel]=coupleRxnList2Rxn(pairedModel,pairedModel.rxns(strmatch(nameTagsModels{2,1},pairedModel.rxns)),strcat(infoFile{modelInd(2,i),1},'_',bioID{2,1}));
    pairedModel=useDiet(pairedModel,dietConstraints{1});
    [ParetoFrontier] = computeParetoOptimality(pairedModel,strcat(infoFile{modelInd(1,i),1},'_',bioID{1,1}),strcat(infoFile{modelInd(2,i),1},'_',bioID{2,1}),'dinc',0.001,'FVAflag',false);
end
%% 
% Can you interpret the shapes of the five Pareto frontiers that were computed? 
% Are there microbe pairs that are always competing with each other? Are there 
% pairs in which one microbe can benefit the other at certain points in the curve 
% and vice versa?



%% so once i have determined growth, what's next

% need some in vitro growth data for parameterizing gLV and determining
% sign of gamma 

% can make my weakly informative prior just the OLS solution of the gLV
% system of ODEs


