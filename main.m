%%
% *Identifying immunologically-vulnerable regions of the HCV E2 
% glycoprotein and broadly neutralizing antibodies that target them*

%% Setting up paths (of functions and data files required) and necessary parameters
clear all;close all;clc

addpath data
addpath functions

run startup.m

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

run_scripts = 0;
%1 = yes, run the scripts. Note that this would take very long time...
%0 = no, use saved data.

%% 1. Preprocessing of the E2 sequence data downloaded from NCBI

% HCV genotype 1a E2 sequences were downloaded from the Los Alamos National 
% Laboratory (LANL) HCV sequence database (https://hcv.lanl.gov; 
% accessed Sep. 25, 2017)
% Header format: Name.PATid.Accession.Subtype.Country.SamplingYear

inputfile = 'HCV_E2_aa_align_sequences_1909.fasta';
[X,w] = preprocessing_data(inputfile);

%% 2. Inferring the model representing the prevalence landscape of E2 using MPF-BML

% Code for running MPF-BML is freely available at <https://github.com/raymondlouie/MPF-BML>. 

%% 
% MPF-BML input: MSA matrix "X" and weight vector "w".

%% 
% MPF-BML output: Maximum entropy model parameters (fields and couplings)

load MPF_BML_output.mat

msa_aa_ex = msa_bin;

[~,L] = size(msa_aa);
% L is the total number of residues in E2
[~,L_ex] = size(msa_aa_ex);
% L_ex is the length of the Potts binary extended MSA
L_mut = length(ind_non_conserve);
% L_mut is the number of non-conserve (mutating) residues in E2

H = reshape(J_mat,L_ex,L_ex); 
% Model-paramters matrix which consists of inferred fields on its diagonal
% and inferred couplings on the upper triangular matrix.

phi_curr = num_mutants_combine_array; 
%Number of mutants at each residue

mutant_order = amino_single_combine_array; 
%Identity of mutants at each residue

phi_cumulative(1) = phi_curr(1);
for kk = 2:length(mutant_order)
    phi_cumulative(kk) = sum(phi_curr(1:kk));
end
%Cumulative sum of the number of mutants at each residue


% Generating samples from the inferred model using a Markov Chain Monte Carlo (MCMC) method

if run_scripts == 1
    samples_MCMC = generate_samples_MCMC(H,msa_aa_ex,phi_curr,phi_cumulative,...
        L_mut,L_ex,w);
else
    load samples_MCMC_E2_99900.mat
end

%% 3. Statistical validation of the inferred model 

% Single mutant probability, double mutant probability, and distribution of
% the number of mutations
model_statistical_validation(msa_aa_ex,w,samples_MCMC,phi_curr,L_mut);

% Triple mutant probability
if run_scripts == 1
    [p3_data, p3_sampler] = compute_triple_mutant_probability(msa_aa,...
        phi_curr,mutant_order,ind_conserve,samples_MCMC,w);
else
    load three_point_correlations_E2_MSA.mat
    load three_point_correlations_E2_model.mat
end
plot_triple_mutant_probability(p3_data,p3_sampler)

%% 4. Biological validation of the inferred fitness landscape

%% 
% In silico predicted energy vs in vitro replicative fitness measurements

% Energies predicted using the inferred model which incorporates the
% effects of mutations at individual residues as well as  interactions
% between mutations at different residues

method = 1;
energy_vs_fitness(method, msa_aa_ex, msa_aa, phi_curr, mutant_order, ...
    H, ind_non_conserve);

% Energy predicted using the model which incorporates only the effects of
% mutations at individual residues

method = 2;
energy_vs_fitness(method, msa_aa_ex, msa_aa, phi_curr, mutant_order, ...
    H, ind_non_conserve);

%% 
% Predicting the easiest compensatory pathway in the H77 background with
% N417S mutation

escape_mutation_N417S_Q444R(msa_aa, phi_curr, mutant_order, ...
    H, ind_non_conserve)

%% 5. Biological validation of the evolutionary "escape time" metric and design of a binary classifier 
% Comparison of escape times associated with known escape mutations and
% the mutations at other residues

% Computing ecscape times (t_e^i) requires to run evolutionay simulation
% model, which is computationally quite complicated
% We ran the code "WF_E2_script_Ng_b_beta.m" for each residue on a 
% supercomputer using 100 nodes, each comprising 24 cores, 
% and calculated the escape time score for each residue.
% We also computed AUC and a escape time cut-off based on the escape time
% scores predicted for the experimentally/clinically observed escape
% mutations.

% First, we load the calculated escape time scores 

load escape_time_E2.mat mean_escape_time

compare_deltaE_known_escape_mutations(mean_escape_time,L)


%% 
% Comparison of escape times associated with mutations at exposed residues
% those at buried ones

compare_deltaE_exposed_buried_residues(mean_escape_time)

%% 
% Comparison of escape times associated with T cell epitopes associated
% with spontaneous cleanrance

compare_deltaE_Tcell_epitopes(mean_escape_time,ind_non_conserve)

%% 6. Mapping of escape times on the available crystal structure of E2

mapping_crystal_structure_pymol(mean_escape_time)
% The above code generates the text file required in Pymol for creating the
% heatmap. The code to run in Pymol is provided in "Pymol_code.docx" file.

%% 7. Comparison of escape times associated with different regions targeted by HmAbs

compare_deltaE_antigenic_domains(mean_escape_time)

%% 8. Comparison of escape times associated with the binding residues of known HmAbs (RB<=20)

compare_deltaE_HmAbs(mean_escape_time)

%% 9. Comparison of escape times associated with the binding residues of known HmAbs (RB<=40)
compare_deltaE_HmAbs_RB40(mean_escape_time)

%% 10. Comparison of escape times associated with the binding residues of known HmAbs (selective alanine scanning)
compare_deltaE_HmAbs_selective(mean_escape_time)
