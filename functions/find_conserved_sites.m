function [indices_conserved_sites,no_of_conserved_sites] = find_conserved_sites(msa)

% Code for identifying 100% conserved sites
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-05

% INPUT
%   - msa: a matrix consisting of E2 sequences in each row
% 
% OUTPUTS
%   - indices_conserved_sites: list of 100% conserved sites
%   - no_of_conserved_sites: number of 100% conserved sites

%%
[ns,~]=size(msa); %ns=no.of seq, ls=length of seq
profile_seq = seqprofile(msa);
[~,pos_wt] = max(profile_seq);
Cseq = int2aa(pos_wt);
Cseq_mtrx = repmat(Cseq,ns,1);
msa_binary = double(msa==Cseq_mtrx); %replacing mutation by 0 and no mutation by 1
freq_single_mutation = mean(msa_binary);
site_freq_no_mutation = find(freq_single_mutation==1);

indices_conserved_sites = site_freq_no_mutation;

no_of_conserved_sites = length(indices_conserved_sites);