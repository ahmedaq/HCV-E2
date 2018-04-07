function dE2 = compute_deltaE(samples_MCMC,msa_aa,phi_curr,mutant_order,...
    ind_non_conserve,ind_conserve,H)

% Code for computing fitness cost (deltaE_i) for each resdiue i in E2
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%
run startup.m

msa_aa_mut = msa_aa(:,ind_non_conserve);

H = convertRayFLParamsToJohnParamsFormat(H);

phi_cumulative(1) = phi_curr(1);
for kk = 2:length(mutant_order)
    phi_cumulative(kk) = sum(phi_curr(1:kk));
end
phi_cum = [0 phi_cumulative]; %for proper indexing in loop

%Generating bin_matrix (definition of binary encoding)
for aba=1:30
    temp_matrix=[];
    temp_matrix = fliplr(eye(aba));
    temp_matrix = [zeros(1,size(temp_matrix,2)) ; temp_matrix ];
    bin_matrix{aba} =temp_matrix;
end

[~,ls_mut] = size(msa_aa_mut);
[~,ls] = size(msa_aa);

ns_MCMC = size(samples_MCMC,1);
%%

parfor kk = 1:ns_MCMC
    EnergySeq_MCMC(kk) = calcSeqEnergy(samples_MCMC(kk,:),H);
end
% save EnergySeq_MCMC_99900 EnergySeq_MCMC

% load EnergySeq_MCMC

probSeq_MCMC = exp(-EnergySeq_MCMC)./sum(exp(-EnergySeq_MCMC));
[a,b] = sort(probSeq_MCMC,'descend');
%%


for kk = 1:ls_mut
    mut_pos = phi_cum(kk)+1:phi_cum(kk+1);
    seqs_MCMC_with_wt_mut = find(sum(samples_MCMC(:,mut_pos),2)==0);
    
    for nn = 1:phi_curr(kk) %loop for mutants at this position
        temp = 0;
        
        parfor mm = 1:length(seqs_MCMC_with_wt_mut) %loop for each seq with wt at this position
            
            SeqWithMut = samples_MCMC(seqs_MCMC_with_wt_mut(mm),:);
            SeqWithMut(phi_cum(kk)+1:phi_cum(kk+1)) = bin_matrix{phi_curr(kk)}(nn+1,:);
            EnergySeqWithMut = calcSeqEnergy(SeqWithMut,H);
            
            temp = temp + (EnergySeq_MCMC(seqs_MCMC_with_wt_mut(mm)) - EnergySeqWithMut) ...
                * probSeq_MCMC(seqs_MCMC_with_wt_mut(mm));
        end
        
        deltaE{kk}(nn) = temp;
    end
    kk
end

% save deltaE_all_99900 deltaE

%% Averaging over all mutants at each residue

dE = zeros(1,ls_mut);
parfor kk = 1:ls_mut
    for nn = 1:phi_curr(kk)
        dE(kk) = dE(kk) + deltaE{kk}(nn)*exp(-deltaE{kk}(nn))/sum(exp(-deltaE{kk}));
    end
end

%% Constructing dE vector for all positions (including the conserved sites)

for kk = 1:ls
    if ismember(kk,ind_conserve)
        %replacing change in energy of conserved sites with the minimum energy difference observed for a mutation
        dE2(kk) = -min(dE); 
    else
        dE2(kk) = -dE(rev_translation_indices(kk,ind_non_conserve));
    end
end

% save dE_99900 dE

% save data_fitnessCosts_E2_99900 