function samples_MCMC = generate_samples_MCMC(H,msa_aa_ex,phi_curr,phi_cumulative,...
    protein_length_aa,protein_length_ex,w)

% Code for generating MCMC samples
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%% Parameters

thin = 1e3;%1e3;
burnin = 1e4;%1e5
nosim = 1e8; %5*1e7; % 1e8;

MCsweepLength = 1e4; % 2*1e4; % 5*1e4; % 1e5;
numParallel = 1; %number of parallel chains
betaArray = 1; %[1 0.85 0.7 0.55]; %[1 0.9 0.8 0.7]; % [1 0.85 0.7 0.55];%1

% this is the seed sequence of the MCMC sampler
seedSeq = rand(1,protein_length_aa) > 0.7;
% seedSeq = zeros(1,protein_length_aa);

param_verifyMCMC_TC_PT = cell(1,12);
param_verifyMCMC_TC_PT{1} = H;
param_verifyMCMC_TC_PT{2} = protein_length_ex;
param_verifyMCMC_TC_PT{3} = msa_aa_ex; % comment out this line
param_verifyMCMC_TC_PT{4} = phi_curr;
param_verifyMCMC_TC_PT{5} = phi_cumulative;
param_verifyMCMC_TC_PT{6} = protein_length_aa;
param_verifyMCMC_TC_PT{9} = [thin burnin nosim];
param_verifyMCMC_TC_PT{10} = [MCsweepLength numParallel betaArray];
param_verifyMCMC_TC_PT{11} = seedSeq;
param_verifyMCMC_TC_PT{12} = w;

%% Sampler

samples_MCMC = verifyMCMC_TC_PT_E2(param_verifyMCMC_TC_PT);
