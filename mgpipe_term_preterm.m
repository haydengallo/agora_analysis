%%% Hayden Gallo
%%% Bucci Lab
%%% 6/17/24 
%%% Running mgPipe from Agora/cobratoolbox


%% Initializing cobratoolbox and setting gurobi solver as default

initCobraToolbox

% gurobi specifiec, but can also use ibm cmplx if installed, gurobi is
% often easier to get working
solverOK=changeCobraSolver('ibm_cplex','all');

%% Changing working directory based on machine being run on, either linux or osx

name ='donatello'

%[ret, name] = system('hostname');

% Preparing input data


if (name == 'donatello')
    % path to .mat files
    modPath = '/lacie_donatello/hgallo/cobra_analysis/AGORA_2/AGORA2_mat/';
    % changing working directory
    cd '/lacie_donatello/hgallo/cobra_analysis/';
    % path to abundance file 
    abundance_file = '/lacie_donatello/hgallo/cobra_analysis/agora_term_preterm_data/AMANHI_P_concat_species_level.csv';
    
    %taxa info table 
    taxTable = '/lacie_donatello/hgallo/cobra_analysis/AGORA2_infoFile_AMANHI_P.xlsx';

else 

    % setting directory where .mat files are stored for agora models

    modPath = '/Users/haydengallo/cobratoolbox/AGORA-2/AGORA_2_mat/';

    % then we want to normalize coverage and set cutoff 

    cd '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/';

    % determine which file to use 

    abundance_file = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/AMANHI_P_concat_species_level.csv';

    % taxa info table
    taxTable = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/AGORA2_infoFile_AMANHI_P.xlsx';
end
%% Read in file with abundances and specified .mat models as index


cutoff = 0.0001;
[normalizedCoverage,abundancePath] = normalizeCoverage(abundance_file);

% if needing to create pan models (i.e. not running at the strain level)
% specify yes here

panmodelsneeded = 'yes';


%% Running panmodels if necessary

% set number of workers for parallel processing 

numWorkers = 12;
%%
% creating pan models if needed

if(panmodelsneeded == 'yes')
    
    % specify folder to store panmodels

    panPath=[pwd filesep 'panSpeciesModels_AMANHI_P'];

    % specify level of pan models

    taxonLevel='Species';

    % function to create panmodels

    createPanModels(modPath,panPath,taxonLevel,numWorkers,taxTable);

end



%% Setting parameters for mgPipe


% set diet to be used for environmental basis of simulation
% diets are located in
% ~/cobratoolbox/papers/2018_microbiomemodelingToolbox/input

% dietFilePath='/Users/haydengallo/cobratoolbox/papers/2018_microbiomemodelingToolbox/input/AverageEuropeanDiet';
dietFilePath ='/home/hgallo/cobratoolbox/papers/2018_microbiomemodelingToolbox/input/AverageEuropeanDiet';

% set file path to normalized abundance table 

% this is old can probably get rid of it 
%abunFilePath='/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/normCoverageReduced.csv';

abunFilePath ='/lacie_donatello/hgallo/cobra_analysis/normalizedCoverage.csv';

% this determines if FVA is computed

computeProfiles = true;

% path to csv file with stratification information 

infoFilePath = '/lacie_donatello/hgallo/cobra_analysis/agora_term_preterm_data/sampInfo_AMANHI_P.csv';



%% Running mgPipe function 


panmodPath = '/lacie_donatello/hgallo/cobra_analysis/panSpeciesModels_AMANHI_P';
% panmodPath = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/panSpeciesModels_AMANHI_P';

[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary] = initMgPipe(panmodPath, abundancePath, computeProfiles, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'numWorkers', numWorkers);


%% Correlation between fluxes and abundances 

% don't need this either should be able to just use taxTable 
%taxInfo = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/AGORA2_infoFile_AMANHI_P.xlsx';

%fluxPath = '/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/Results/Results/inputDiet_net_secretion_fluxes.csv';

fluxPath = '/lacie_donatello/hgallo/cobra_analysis/Results/inputDiet_net_secretion_fluxes.csv';

%abunFilePath ='/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/normalizedCoverage_for_corr_only.csv'

abunFilePath = '/lacie_donatello/hgallo/cobra_analysis/normalizedCoverage_corr_only.csv';

corrMethod = 'Spearman';
[FluxCorrelations, PValues, TaxonomyInfo] = correlateFluxWithTaxonAbundance(abunFilePath, fluxPath, taxTable, corrMethod);

%% fluxes against organism abundances 
plotFluxesAgainstOrganismAbundances(abunFilePath,fluxPath);

%% Stratification of samples


[init, netSecretionFluxes, netUptakeFluxes, Y, modelStats, summary, statistics] = initMgPipe(modPath, abunFilePath, computeProfiles, 'dietFilePath', dietFilePath, 'infoFilePath', infoFilePath, 'numWorkers', numWorkers);

%% Analysis of mgPipe output, statistical tests to determine if fluxes or uptakes are different between groups (term vs. pre-term)

% Plotting of fluxes 


% metadata file where grouping is specified



%infoFilePath='/Users/haydengallo/Documents/Bucci_Lab/cobratoolbox_analysis/agora_term_preterm_data/sampInfo_AMANHI_P.csv';


sampleGroupHeaders={'Group'};


resPath = '/lacie_donatello/hgallo/cobra_analysis/Results';


analyzeMgPipeResults(infoFilePath,resPath, 'sampleGroupHeaders', sampleGroupHeaders);


%% Determining strain-level contributions to fluxes and metabolite flows 

% location of models constrained by diet 

constrModPath = [resPath filesep 'Diet'];

% list of specified metabolites to analyze, if none listed, then all
% exchanged metabolites will be analyzed 

%metList = {'ac','for','succ', 'but'};

% calling actual function to determine strain contributions 

[minFluxes,maxFluxes,fluxSpans] = predictMicrobeContributions(constrModPath,'numWorkers', numWorkers);


%%

% path where the conributions are located

contPath = '/lacie_donatello/hgallo/cobra_analysis/Contributions';


% running actual analysis 

analyzeMgPipeResults(infoFilePath,contPath, 'sampleGroupHeaders', sampleGroupHeaders);











