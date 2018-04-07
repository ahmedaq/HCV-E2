
function samples_MCMC_double = verifyMCMC_TC_PT_E2(param_verifyMCMC_TC_PT)

% Code for running MCMC sampler
% 
% Written by: Raymond Louie and Muhammad Saqib Sohail
% Edited by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%
run startup.m

J_MINFLOW_mat = param_verifyMCMC_TC_PT{1};
total_length = param_verifyMCMC_TC_PT{2};
msa_aa_ex = param_verifyMCMC_TC_PT{3};
phi_curr = param_verifyMCMC_TC_PT{4};
phi_cumulative = param_verifyMCMC_TC_PT{5};
protein_length_aa = param_verifyMCMC_TC_PT{6};
temp12 = param_verifyMCMC_TC_PT{9};
temp22 = param_verifyMCMC_TC_PT{10};
seedSeq = param_verifyMCMC_TC_PT{11};
weight = param_verifyMCMC_TC_PT{12};

thin = temp12(1);
burnin = temp12(2);
nosim = temp12(3);

MCsweepLength = temp22(1);
numParallel = temp22(2);
betaArray = temp22(3:end);


clear param_verifyMCMC;


%% verfication

random_array_sweep = rand(1,nosim/MCsweepLength*(numParallel - 1));
number_samples = ceil((nosim-burnin)/thin);


J_MINFLOW_mat_array = J_MINFLOW_mat(:);
t_samp = tic();

random_site_array =  randi([1 protein_length_aa],1,nosim*numParallel); % choose the random site

num_amino_array = phi_curr(random_site_array);
unique_amino = unique(num_amino_array);
for cbc=1:length(unique_amino)
    num_amino = unique_amino(cbc);
    ind = find(num_amino_array==unique_amino(cbc));
    rand_amino_array(ind) = randi([1 num_amino],1,length(ind));
end
random_array = rand(1,nosim*numParallel);
t_samp = toc(t_samp);
fprintf( 'Random generation in %f seconds \n', t_samp );

% Run MCMC

totalnosample=0;
t_samp = tic();
samples_MCMC_double = zeros(number_samples,total_length);

aba=1;

curr_vector = seedSeq;

mex PT_final_chain_1.c

% C code written by Muhammad Saqib Sohail
[doublemutant,nosample,energyAll,numMutAll,Vall,samples_MCMC ]= PT_final_chain_1(random_array,...
    random_site_array,rand_amino_array,curr_vector,J_MINFLOW_mat_array,total_length,...
    nosim,phi_cumulative,phi_curr,burnin,thin,curr_vector,number_samples, numParallel, ...
    betaArray(1:numParallel),MCsweepLength, random_array_sweep);

samples_MCMC_double((aba-1)*number_samples+1:(aba-1)*number_samples+number_samples,:)=...
    reshape(samples_MCMC,total_length,nosample)';
%energyVecs = reshape(energyAll, nosim, numParallel); % eaach column is energy of each chain
%numMutVecs = reshape(numMutAll, nosim, numParallel); % eaach column is magnatization of each chain
