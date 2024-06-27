%%% Hayden Gallo
%%% Bucci Lab
%%% 6/17/24 
%%% Running mgPipe from Agora/cobratoolbox


%% Initializing cobratoolbox and setting gurobi solver as default

initCobraToolbox

% gurobi specifiec, but can also use ibm cmplx if installed, gurobi is
% often easier to get working
solverOK=changeCobraSolver('gurobi','LP');


%% Preparing input data for mgPipe

% setting directory where .mat files are stored for agora models

modPath = '/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/';



%% Read in file with abundances and specified .mat models as index
% then we want to normalize coverage and set cutoff 

cd '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/';

% determine which file to use 

abundance_file = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/AMANHI_P_concat_species_level.csv';


cutoff = 0.0001;
[normalizedCoverage,abundancePath] = normalizeCoverage(abundance_file);

% if needing to create pan models (i.e. not running at the strain level)
% specify yes here

panmodelsneeded = 'yes';


%% Running panmodels if necessary

% set number of workers for parallel processing 

numWorkers = 8;
%%
% creating pan models if needed

if(panmodelsneeded == 'yes')
    
    % specify folder to store panmodels

    panPath=[pwd filesep 'panSpeciesModels_AMANHI_P'];

    % specify level of pan models

    taxonLevel='Species';

    % function to create panmodels

    createPanModels(modPath,panPath,taxonLevel,'taxTable','/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/AGORA2_infoFile_AMANHI_P.xlsx');

end



%% Setting parameters for mgPipe


% set diet to be used for environmental basis of simulation
% diets are located in
% ~/cobratoolbox/papers/2018_microbiomemodelingToolbox/input

dietFilePath='/Users/haydengallo/cobratoolbox/papers/2018_microbiomemodelingToolbox/input/AverageEuropeanDiet';


% set file path to normalized abundance table 

abunFilePath='/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/normCoverageReduced.csv';

% this determines if FVA is computed

computeProfiles = true;

% path to csv file with stratification information 

infoFilePath = '';



%% Running mgPipe function 


panmodPath = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/panSpeciesModels_AMANHI_P';

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] = initMgPipe(panmodPath, abundancePath, computeProfiles, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'numWorkers', numWorkers);


%% Correlation between fluxes and abundances 

taxInfo = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/AGORA2_infoFile_AMANHI_P.xlsx';

fluxPath = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/Results/Results/inputDiet_net_secretion_fluxes.csv';

abunFilePath ='/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/normalizedCoverage_for_corr_only.csv'

corrMethod = 'Spearman';
[FluxCorrelations, PValues, TaxonomyInfo] = correlateFluxWithTaxonAbundance(abunFilePath, fluxPath, taxInfo, corrMethod);

%% fluxes against organism abundances 
plotFluxesAgainstOrganismAbundances(abunFilePath,fluxPath,metList);

%% Stratification of samples


[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary, statistics] = initMgPipe(modPath, abunFilePath, computeProfiles, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'numWorkers', numWorkers);

%% Analysis of mgPipe output, statistical tests to determine if fluxes or uptakes are different between groups (term vs. pre-term)

% Plotting of fluxes 


% metadata file where grouping is specified



infoFilePath='/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/sampInfo_AMANHI_P.csv';


sampleGroupHeaders={'Group'};


resPath = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/Results/Results';


analyzeMgPipeResults(infoFilePath,resPath, 'sampleGroupHeaders', sampleGroupHeaders);


%% Determining strain-level contributions to fluxes and metabolite flows 

% location of models constrained by diet 

constrModPath = [resPath filesep 'Diet'];

% list of specified metabolites to analyze, if none listed, then all
% exchanged metabolites will be analyzed 

metList = {'ac','for','succ', 'but'};

% calling actual function to determine strain contributions 

[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(constrModPath,'metList', metList,'numWorkers', numWorkers);

% path where the conributions are located

contPath = [pwd filesep 'Contributions'];


% running actual analysis 

%analyzeMgPipeResults(infoFilePath,contPath, 'sampleGroupHeaders', sampleGroupHeaders);











