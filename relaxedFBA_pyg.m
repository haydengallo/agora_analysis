%% Hayden Gallo
% Bucci Lab
% Cobratoolbox Modeling
% 7/29/24
% Determining reduced cost of different medias to obtain growth from F.
% prau and P. copri


%% initialize cobratoolbox

initCobraToolbox


%% Load different models 

p_copri_model_alt = readCbModel('/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/Prevotella_copri_CB7_DSM_18205.mat');
f_prau_model = readCbModel('/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/Faecalibacterium_prausnitzii_ERR1022327.mat');
combined_model = readCbModel('/Users/haydengallo/cobratoolbox/Results_p_copri_f_prau/pairedModel_Prevotella_copri_ERR1022397_Faecalibacterium_prausnitzii_ERR1022484.mat');


% load different models of F. prau and P. copri that constrains need to be
% examined

[~,infoFile,~]=xlsread('/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/p_copri_f_prau_only.xlsx');

% grab the model names so that can loop through and constrain later on 
modelList=infoFile(2:end,1);
modPath = '/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/';

%% Diets/Medias

%pyg_media = {'EX_mobd(e)','-3.07E-04','1000';'EX_ac(e)','-0.0','1000';'EX_n2(e)','-227.8571429','1000';'EX_na1(e)','-20.51163798','1000';'EX_cl(e)','-5.941065976','1000';'EX_ca2(e)','-0.173941043','1000';'EX_fe2(e)','-0.016053362','1000';'EX_mg2(e)','-0.474477191','1000';'EX_k(e)','-35.39748582','1000';'EX_so4(e)','-0.710983412','1000';'EX_pi(e)','-18.29826648','1000';'EX_ala_L(e)','-16.89227108','1000';'EX_arg_L(e)','-5.338568575','1000';'EX_asn_L(e)','-1.286718791','1000';'EX_asp_L(e)','-10.81893313','1000';'EX_Lcystin(e)','-0.187272152','1000';'EX_glu_L(e)','-18.56754922','1000';'EX_gln_L(e)','-0.205274178','1000';'EX_gly(e)','-20.3151851','1000';'EX_his_L(e)','-3.319218598','1000';'EX_ile_L(e)','-8.195159139','1000';'EX_leu_L(e)','-11.2826377','1000';'EX_lys_l(e)','-9.884397018','1000';'EX_met_L(e)','-2.144657123','1000';'EX_phe_L(e)','-6.083829725','1000';'EX_pro_L(e)','-11.89938505','1000';'EX_ser_L(e)','-4.424652451','1000';'EX_thr_L(e)','-3.567830759','1000';'EX_trp_L(e)','-0.685504997','1000';'EX_tyr_L(e)','-1.683306566','1000';'EX_val_L(e)','-10.37149589','1000';'EX_glc_D(e)','-27.7537245498346','1000';'EX_hco3(e)','-1.190379826','1000';'EX_phyQ(e)','-0.002172408','1000';'EX_etoh(e)','-3.237913066','1000';'EX_pheme(e)','-0.0076693','1000';'EX_oh1(e)','-0.099987502','1000';'EX_cys_L(e)','-2.846894039','1000';'EX_M02144(e)','-2.846894039','1000';'EX_h2o(e)','-55013.2623','1000';'EX_h(e)','-6.30957E-05','1000'};
pyg_test = {'EX_n2(e)','-227.8571429','1000';'EX_na1(e)','-20.51163798','1000';'EX_cl(e)','-5.941065976','1000';'EX_ca2(e)','-0.173941043','1000';'EX_fe2(e)','-0.016053362','1000';'EX_mg2(e)','-0.474477191','1000';'EX_k(e)','-35.39748582','1000';'EX_so4(e)','-0.710983412','1000';'EX_pi(e)','-18.29826648','1000';'EX_ala_L(e)','-16.89227108','1000';'EX_arg_L(e)','-5.338568575','1000';'EX_asn_L(e)','-1.286718791','1000';'EX_asp_L(e)','-10.81893313','1000';'EX_Lcystin(e)','-0.187272152','1000';'EX_glu_L(e)','-18.56754922','1000';'EX_gln_L(e)','-0.205274178','1000';'EX_gly(e)','-20.3151851','1000';'EX_his_L(e)','-3.319218598','1000';'EX_ile_L(e)','-8.195159139','1000';'EX_leu_L(e)','-11.2826377','1000';'EX_lys_l(e)','-9.884397018','1000';'EX_met_L(e)','-2.144657123','1000';'EX_phe_L(e)','-6.083829725','1000';'EX_pro_L(e)','-11.89938505','1000';'EX_ser_L(e)','-4.424652451','1000';'EX_thr_L(e)','-3.567830759','1000';'EX_trp_L(e)','-0.685504997','1000';'EX_tyr_L(e)','-1.683306566','1000';'EX_val_L(e)','-10.37149589','1000';'EX_glc_D(e)','-27.7537245498346','1000';'EX_hco3(e)','-1.190379826','1000';'EX_phyQ(e)','-0.002172408','1000';'EX_etoh(e)','-3.237913066','1000';'EX_pheme(e)','-0.0076693','1000';'EX_oh1(e)','-0.099987502','1000';'EX_cys_L(e)','-2.846894039','1000';'EX_M02144(e)','-2.846894039','1000';'EX_h2o(e)','-55013.2623','1000';'EX_h(e)','-6.30957E-05','1000'};

pyg_test_metals_add = {'EX_zn2(e)','-1','1000';'EX_mn2(e)','-1','1000';'EX_fe3(e)','-1','1000';'EX_cu2(e)','-1','1000';'EX_cobalt2(e)','-1','1000';'EX_n2(e)','-227.8571429','1000';'EX_na1(e)','-20.51163798','1000';'EX_cl(e)','-5.941065976','1000';'EX_ca2(e)','-0.173941043','1000';'EX_fe2(e)','-0.016053362','1000';'EX_mg2(e)','-0.474477191','1000';'EX_k(e)','-35.39748582','1000';'EX_so4(e)','-0.710983412','1000';'EX_pi(e)','-18.29826648','1000';'EX_ala_L(e)','-16.89227108','1000';'EX_arg_L(e)','-5.338568575','1000';'EX_asn_L(e)','-1.286718791','1000';'EX_asp_L(e)','-10.81893313','1000';'EX_Lcystin(e)','-0.187272152','1000';'EX_glu_L(e)','-18.56754922','1000';'EX_gln_L(e)','-0.205274178','1000';'EX_gly(e)','-20.3151851','1000';'EX_his_L(e)','-3.319218598','1000';'EX_ile_L(e)','-8.195159139','1000';'EX_leu_L(e)','-11.2826377','1000';'EX_lys_l(e)','-9.884397018','1000';'EX_met_L(e)','-2.144657123','1000';'EX_phe_L(e)','-6.083829725','1000';'EX_pro_L(e)','-11.89938505','1000';'EX_ser_L(e)','-4.424652451','1000';'EX_thr_L(e)','-3.567830759','1000';'EX_trp_L(e)','-0.685504997','1000';'EX_tyr_L(e)','-1.683306566','1000';'EX_val_L(e)','-10.37149589','1000';'EX_glc_D(e)','-27.7537245498346','1000';'EX_hco3(e)','-1.190379826','1000';'EX_phyQ(e)','-0.002172408','1000';'EX_etoh(e)','-3.237913066','1000';'EX_pheme(e)','-0.0076693','1000';'EX_oh1(e)','-0.099987502','1000';'EX_cys_L(e)','-2.846894039','1000';'EX_M02144(e)','-2.846894039','1000';'EX_h2o(e)','-55013.2623','1000';'EX_h(e)','-6.30957E-05','1000'};
pyg_test_all_add = {'EX_q8(e)','-1','1000';'EX_pydx(e)','-1','1000';'EX_nac(e)','-1','1000';'EX_lys_L(e)','-1','1000';'EX_hxan(e)','-1','1000';'EX_ade(e)','-1','1000';'EX_thymd(e)','-1','1000';'EX_thm(e)','-1','1000';'EX_ribflv(e)','-1','1000';'EX_pnto_R(e)','-1','1000';'EX_nac(e)','-1','1000';'EX_fol(e)','-1','1000';'EX_zn2(e)','-1','1000';'EX_mn2(e)','-1','1000';'EX_fe3(e)','-1','1000';'EX_cu2(e)','-1','1000';'EX_cobalt2(e)','-1','1000';'EX_n2(e)','-227.8571429','1000';'EX_na1(e)','-20.51163798','1000';'EX_cl(e)','-5.941065976','1000';'EX_ca2(e)','-0.173941043','1000';'EX_fe2(e)','-0.016053362','1000';'EX_mg2(e)','-0.474477191','1000';'EX_k(e)','-35.39748582','1000';'EX_so4(e)','-0.710983412','1000';'EX_pi(e)','-18.29826648','1000';'EX_ala_L(e)','-16.89227108','1000';'EX_arg_L(e)','-5.338568575','1000';'EX_asn_L(e)','-1.286718791','1000';'EX_asp_L(e)','-10.81893313','1000';'EX_Lcystin(e)','-0.187272152','1000';'EX_glu_L(e)','-18.56754922','1000';'EX_gln_L(e)','-0.205274178','1000';'EX_gly(e)','-20.3151851','1000';'EX_his_L(e)','-3.319218598','1000';'EX_ile_L(e)','-8.195159139','1000';'EX_leu_L(e)','-11.2826377','1000';'EX_lys_l(e)','-9.884397018','1000';'EX_met_L(e)','-2.144657123','1000';'EX_phe_L(e)','-6.083829725','1000';'EX_pro_L(e)','-11.89938505','1000';'EX_ser_L(e)','-4.424652451','1000';'EX_thr_L(e)','-3.567830759','1000';'EX_trp_L(e)','-0.685504997','1000';'EX_tyr_L(e)','-1.683306566','1000';'EX_val_L(e)','-10.37149589','1000';'EX_glc_D(e)','-27.7537245498346','1000';'EX_hco3(e)','-1.190379826','1000';'EX_phyQ(e)','-0.002172408','1000';'EX_etoh(e)','-3.237913066','1000';'EX_pheme(e)','-0.0076693','1000';'EX_oh1(e)','-0.099987502','1000';'EX_cys_L(e)','-2.846894039','1000';'EX_M02144(e)','-2.846894039','1000';'EX_h2o(e)','-55013.2623','1000';'EX_h(e)','-6.30957E-05','1000'};



%% Apply diet to P. copri 


p_copri_alt_pyg = useDiet(p_copri_model_alt, pyg_test);

%%
excludedReactions = ~startsWith(p_copri_model_alt.rxns, 'EX_');


%%

param = struct();
%param.internalRelax = 0;  % Allow to relax bounds on all internal reactions
%param.exchangeRelax = 2;  % Do not allow to relax bounds on exchange reactions
%param.steadyStateRelax = 0;
param.excludedReactions = excludedReactions;

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
    %if ~isempty(rxnIndex)
    %    excludedReactionLB(rxnIndex) = true; % Set the corresponding index to true
    %end


%%
% Set all lower bounds to zero
%p_copri_alt_pyg.lb(:) = 0;


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

p_copri_alt_pyg.lb


%%
%p_copri_alt_pyg = useDiet(p_copri_model_alt, pyg_test);
p_copri_alt_pyg = changeRxnBounds(p_copri_alt_pyg, 'biomass525', 0.001, 'l');

%%

% Add the excludedReactionLB to the param structure
%param.excludedReactionLB = excludedReactionLB;


%%



[solution, relaxedModel] = relaxedFBA(p_copri_alt_pyg, param);

% Check if a solution was found
if isempty(solution.v)
    disp('relaxedFBA could not find a solution');
else
    disp('relaxedFBA found a solution');
    
    % Apply relaxed constraints
    relaxedModel.lb = relaxedModel.lb - solution.p;
    relaxedModel.ub = relaxedModel.ub + solution.q;
    
    % Try to optimize with relaxed constraints
    relaxedSolution = optimizeCbModel(relaxedModel);
    
    if isempty(relaxedSolution.x)
        disp('Model is still infeasible after relaxation');
    else
        disp(['Growth rate after relaxation: ', num2str(relaxedSolution.f)]);
    end
end


%%

disp('Excluded Reactions (expected true for non-EX_ reactions):');
disp(param.excludedReactions);


%%
% Get the bounds from the relaxed model
relaxedLB = relaxedModel.lb;
relaxedUB = relaxedModel.ub;

temp_model_pygLB = p_copri_alt_pyg.lb;
temp_model_pygUB = p_copri_alt_pyg.ub;



% Compare and print the differences
disp('Reactions with changed bounds:');
disp('Reaction ID | P_copri_PYG LB | Relaxed LB | P_copri_PYG UB | Relaxed UB');
for i = 1:length(p_copri_alt_pyg.rxns)
    if temp_model_pygLB(i) ~= relaxedLB(i) || temp_model_pygUB(i) ~= relaxedUB(i)
        fprintf('%s | %.6f | %.12f | %.4f | %.4f\n', ...
            p_copri_alt_pyg.rxns{i}, temp_model_pygLB(i), relaxedLB(i), temp_model_pygUB(i), relaxedUB(i));
    end
end

% Calculate and print statistics
numChanged = sum((temp_model_pygLB ~= relaxedLB) | (temp_model_pygUB ~= relaxedUB));
percentChanged = (numChanged / length(p_copri_alt_pyg.rxns)) * 100;

disp(['Number of reactions with changed bounds: ', num2str(numChanged)]);
disp(['Percentage of reactions with changed bounds: ', num2str(percentChanged), '%']);

% Identify the reactions with the largest relaxations
%[~, indicesLB] = sort(relaxedLB - p_copri_alt_pygLB, 'descend');
%[~, indicesUB] = sort(relaxedUB - p_copri_alt_pygUB, 'descend');

%disp('Top 10 reactions with largest lower bound relaxations:');
%for i = 1:10
%    fprintf('%s: %.4f\n', p_copri_alt_pyg.rxns{indicesLB(i)}, relaxedLB(indicesLB(i)) - p_copri_alt_pygLB(indicesLB(i)));
%end
%
%disp('Top 10 reactions with largest upper bound relaxations:');
%for i = 1:10
%    fprintf('%s: %.4f\n', p_copri_alt_pyg.rxns{indicesUB(i)}, relaxedUB(indicesUB(i)) - p_copri_alt_pygUB(indicesUB(i)));
%end

%% Now testing growth of P copri on media with addition of required metals


p_copri_alt_pyg_add = useDiet(p_copri_model_alt, pyg_test_metals_add);
p_copri_alt_pyg_add = changeRxnBounds(p_copri_alt_pyg_add, 'biomass525', 0.001, 'l');


[test_sol] = optimizeCbModel(p_copri_alt_pyg_add);

%% Ok now need to perform the optimization for F. prau because it doesn't grow on the media with the added metals

f_prau_model_pyg = useDiet(f_prau_model, pyg_test);
f_prau_model_pyg = changeRxnBounds(f_prau_model_pyg, 'bio1', 0.001, 'l');
%%
%excludedReactions = ~startsWith(p_copri_model_alt.rxns, 'EX_');


%%

excludedReactions = ~startsWith(f_prau_model.rxns, 'EX_');

param = struct();
%param.internalRelax = 0;  % Allow to relax bounds on all internal reactions
%param.exchangeRelax = 2;  % Do not allow to relax bounds on exchange reactions
%param.steadyStateRelax = 0;
param.excludedReactions = excludedReactions;






[solution, relaxedModel] = relaxedFBA(f_prau_model_pyg, param);

% Check if a solution was found
if isempty(solution.v)
    disp('relaxedFBA could not find a solution');
else
    disp('relaxedFBA found a solution');
    
    % Apply relaxed constraints
    relaxedModel.lb = relaxedModel.lb - solution.p;
    relaxedModel.ub = relaxedModel.ub + solution.q;
    
    % Try to optimize with relaxed constraints
    relaxedSolution = optimizeCbModel(relaxedModel);
    
    if isempty(relaxedSolution.x)
        disp('Model is still infeasible after relaxation');
    else
        disp(['Growth rate after relaxation: ', num2str(relaxedSolution.f)]);
    end
end






%% F prau
% Get the bounds from the relaxed model
relaxedLB = relaxedModel.lb;
relaxedUB = relaxedModel.ub;

temp_model_pygLB = f_prau_model_pyg.lb;
temp_model_pygUB = f_prau_model_pyg.ub;



% Compare and print the differences
disp('Reactions with changed bounds:');
disp('Reaction ID | P_copri_PYG LB | Relaxed LB | P_copri_PYG UB | Relaxed UB');
for i = 1:length(f_prau_model_pyg.rxns)
    if temp_model_pygLB(i) ~= relaxedLB(i) || temp_model_pygUB(i) ~= relaxedUB(i)
        fprintf('%s | %.6f | %.12f | %.4f | %.4f\n', ...
           f_prau_model_pyg.rxns{i}, temp_model_pygLB(i), relaxedLB(i), temp_model_pygUB(i), relaxedUB(i));
    end
end

% Calculate and print statistics
numChanged = sum((temp_model_pygLB ~= relaxedLB) | (temp_model_pygUB ~= relaxedUB));
percentChanged = (numChanged / length(p_copri_alt_pyg.rxns)) * 100;

disp(['Number of reactions with changed bounds: ', num2str(numChanged)]);
disp(['Percentage of reactions with changed bounds: ', num2str(percentChanged), '%']);

% Identify the reactions with the largest relaxations
%[~, indicesLB] = sort(relaxedLB - p_copri_alt_pygLB, 'descend');
%[~, indicesUB] = sort(relaxedUB - p_copri_alt_pygUB, 'descend');

%disp('Top 10 reactions with largest lower bound relaxations:');
%for i = 1:10
%    fprintf('%s: %.4f\n', p_copri_alt_pyg.rxns{indicesLB(i)}, relaxedLB(indicesLB(i)) - p_copri_alt_pygLB(indicesLB(i)));
%end
%
%disp('Top 10 reactions with largest upper bound relaxations:');
%for i = 1:10
%    fprintf('%s: %.4f\n', p_copri_alt_pyg.rxns{indicesUB(i)}, relaxedUB(indicesUB(i)) - p_copri_alt_pygUB(indicesUB(i)));
%end





%% Now testing growth of F. prau on media with addition of required metals

f_prau_model_pyg_add = useDiet(f_prau_model, pyg_test_all_add);
f_prau_model_pyg_add = changeRxnBounds(f_prau_model_pyg_add, 'bio1', 0.001, 'l');


[test_sol_f_prau] = optimizeCbModel(f_prau_model_pyg_add);


%% simulating P. copri on pyg with all additions 

p_copri_pyg_add = useDiet(p_copri_model_alt, pyg_test_all_add);
p_copri_pyg_add = changeRxnBounds(p_copri_pyg_add, 'biomass525', 0.001, 'l');


[test_sol_p_copri] = optimizeCbModel(p_copri_pyg_add);



%% Here loop through all of the different F. prau and P. copri strains to determine which Exchange reactions are limited 

rxns_to_add ={};


for i = 1:length(modelList)
    
    % load each model

    temp_model = readCbModel([modPath filesep modelList{i} '.mat']);
    
    % constrain each temp model with base pyg media
    temp_model_pyg = useDiet(temp_model, pyg_test);
    
    % Find correct biomass function to constrain
    biomass_rxn = temp_model.rxns{find(strncmp(temp_model.rxns, 'bio', 3)),1};
    disp(biomass_rxn)
 

    temp_model_pyg = changeRxnBounds(temp_model_pyg, biomass_rxn, 0.001, 'l');


    % setting up relaxedFBA reactions to exclude 
    excludedReactions = ~startsWith(temp_model.rxns, 'EX_');
    param = struct();
    param.excludedReactions = excludedReactions;
    
    
    [solution, relaxedModel] = relaxedFBA(temp_model_pyg, param);
    
    % Check if a solution was found
    if isempty(solution.v)
        disp('relaxedFBA could not find a solution');
    else
        disp('relaxedFBA found a solution');
        
        % Apply relaxed constraints
        relaxedModel.lb = relaxedModel.lb - solution.p;
        relaxedModel.ub = relaxedModel.ub + solution.q;
        
        % Try to optimize with relaxed constraints
        relaxedSolution = optimizeCbModel(relaxedModel);
        
        if isempty(relaxedSolution.x)
            disp('Model is still infeasible after relaxation');
        else
            disp(['Growth rate after relaxation: ', num2str(relaxedSolution.f)]);
        end
    end

    % Get the bounds from the relaxed model
    relaxedLB = relaxedModel.lb;
    relaxedUB = relaxedModel.ub;
    
    temp_model_pygLB = temp_model_pyg.lb;
    temp_model_pygUB = temp_model_pyg.ub;
    
    
    
    % Compare and print the differences
    disp('Reactions with changed bounds:');
    disp('Reaction ID | P_copri_PYG LB | Relaxed LB | P_copri_PYG UB | Relaxed UB');
    for i = 1:length(temp_model_pyg.rxns)
        if temp_model_pygLB(i) ~= relaxedLB(i) || temp_model_pygUB(i) ~= relaxedUB(i)
            %rxns_to_add = [rxns_to_add, temp_model_pyg.rxns{i}];
            fprintf('%s | %.6f | %.12f | %.4f | %.4f\n', ...
               temp_model_pyg.rxns{i}, temp_model_pygLB(i), relaxedLB(i), temp_model_pygUB(i), relaxedUB(i));
            rxns_to_add = [rxns_to_add, {temp_model_pyg.rxns{i}}];
        end
    end
    
    % Calculate and print statistics
    numChanged = sum((temp_model_pygLB ~= relaxedLB) | (temp_model_pygUB ~= relaxedUB));
    percentChanged = (numChanged / length(temp_model_pyg.rxns)) * 100;
    
    disp(['Number of reactions with changed bounds: ', num2str(numChanged)]);
    disp(['Percentage of reactions with changed bounds: ', num2str(percentChanged), '%']);



end


disp(unique(rxns_to_add))



%% Loop through all models and make sure there is minimal growth on pyg media with additions of minimal fluxes

for i = 1:length(modelList)
    
    % load each model

    temp_model = readCbModel([modPath filesep modelList{i} '.mat']);
    
    % constrain each temp model with base pyg media
    temp_model_pyg_all = useDiet(temp_model, pyg_test_all_add);
    
    % Find correct biomass function to constrain
    biomass_rxn = temp_model.rxns{find(strncmp(temp_model.rxns, 'bio', 3)),1};
 

    temp_model_pyg_all = changeRxnBounds(temp_model_pyg_all, biomass_rxn, 0.001, 'l');

    temp_model_sol = optimizeCbModel(temp_model_pyg_all)
  
end













