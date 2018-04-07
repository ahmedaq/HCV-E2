function [p3_data, p3_sampler] = ...
    compute_triple_mutant_probability(msa_aa,phi_curr,mutant_order,ind_conserve,samples_MCMC,weight_id)

% Code for computing the triple mutant probabilities in the MSA and the 
% samples generated using the inferred model
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%
msa_aa_mut = msa_aa;
msa_aa_mut(:,ind_conserve) = [];

[ns,ls] = size(msa_aa_mut)

%making the modified msa according to mutant_order

msa_new_data = construct_msa_aa_after_entropy_compression(msa_aa_mut,mutant_order);
%%
%MCMC samples
% load samples_MCMC_E2_99990.mat

phi_cumulative(1) = phi_curr(1);
for kk = 2:length(mutant_order)
    phi_cumulative(kk) = sum(phi_curr(1:kk));
end
phi_cum = [0 phi_cumulative];

for kk = 1:length(phi_curr)
    a = log2(bi2de(fliplr(samples_MCMC(:,phi_cum(kk)+1:phi_cum(kk+1)))))+2;
    a(a==-inf) = 1;
    msa_MCMC_amino(:,kk) = mutant_order{kk}(a);
end

%% calculating 3 point correlation in data
ls = 24;
p3_data = [];
p3_sampler = [];
tic
for kk = 1:ls
    for mm = kk+1:ls
        for nn = mm+1:ls
            p33_data = [];
            p33_sampler = [];
            for kkk = 1:(length(mutant_order{kk})-1)
                for mmm = 1:(length(mutant_order{mm})-1)
                    for nnn = 1:(length(mutant_order{nn})-1)
                        temp_data = msa_new_data(:,[kk mm nn]) == ...
                            repmat([mutant_order{kk}(kkk+1) mutant_order{mm}(mmm+1) ...
                            mutant_order{nn}(nnn+1)],ns,1);
                        p33_data = [p33_data sum((sqrt(weight_id).*temp_data(:,1)) .* ... %important to take care of weight here
                            (sqrt(weight_id).*temp_data(:,2)) .* ...
                            (sqrt(weight_id).*temp_data(:,3)))/sum(weight_id)];
                        
                        temp = msa_MCMC_amino(:,[kk mm nn]) == ...
                            repmat([mutant_order{kk}(kkk+1) mutant_order{mm}(mmm+1) ...
                            mutant_order{nn}(nnn+1)],size(msa_MCMC_amino,1),1);
                        p33_sampler = [p33_sampler sum((temp(:,1).*temp(:,2).*temp(:,3)))/size(msa_MCMC_amino,1)];
                    end
                end
            end
            p3_data = [p3_data fliplr(p33_data)];
            p3_sampler = [p3_sampler fliplr(p33_sampler)];
        end
    end
end
toc