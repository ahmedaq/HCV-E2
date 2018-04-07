function [X,w,indices_conserved_sites] = preprocessing_data(inputfile)

% Code for processing the downloaded sequence data
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-05

% INPUT
%   - inputfile: fasta file of amino acid sequence data of E2 downloaded from LANL
% 
% OUTPUTS
%   - X: a matrix consisting of E2 sequences in each row
%   - w: a vector consisting of patient weight associated with each
%   sequence
%   - indices_conserved_sites: list of 100% conserved sites

%%

run startup.m

%% Processing sequences and associated headers
[header,seqs_aa] = fastaread(inputfile);

M = length(header);

%Constructing matrix MSA

for kk = 1:M
    msa_aa(kk,:) = seqs_aa{kk};
end

[~,L] = size(msa_aa);

% Getting info from header
seqs_noPatID = [];
for kk = 1:M
    
    indx_dot = find(header{kk}=='.');
    if length(indx_dot)==8
        acc_nos{kk} = header{kk}(indx_dot(2)+1:indx_dot(3)-1);
        subtype{kk} = header{kk}(indx_dot(3)+1:indx_dot(4)-1);
        country{kk} = header{kk}(indx_dot(4)+1:indx_dot(5)-1);
        year(kk) = str2double(header{kk}(indx_dot(5)+1:indx_dot(6)-1));
    elseif length(indx_dot)==9
        acc_nos{kk} = header{kk}(indx_dot(3)+1:indx_dot(4)-1);
        subtype{kk} = header{kk}(indx_dot(4)+1:indx_dot(5)-1);
        country{kk} = header{kk}(indx_dot(5)+1:indx_dot(6)-1);
        year(kk) = str2double(header{kk}(indx_dot(6)+1:indx_dot(7)-1));
    
    else 
        acc_nos{kk} = header{kk}(indx_dot(4)+1:indx_dot(5)-1);
        subtype{kk} = header{kk}(indx_dot(5)+1:indx_dot(6)-1);
        country{kk} = header{kk}(indx_dot(6)+1:indx_dot(7)-1);
        year(kk) = str2double(header{kk}(indx_dot(7)+1:indx_dot(8)-1));
    end
    
    patientIDs{kk} = header{kk}(indx_dot(1)+1:indx_dot(2)-1);  
end

%% Finding HXB2 ref seq

for kk = 1:M
    if strcmp(acc_nos{kk},'NC_004102') %HXB2 ref seq
        indx_hxb2_seq = kk;
    end 
end

%% Checking number of amino acids at each position

for kk = 1:L
    no_aa_site(kk) = length(unique(msa_aa(:,kk)));
end 

figure;
bar(1:L,no_aa_site,0.8,'FaceColor',blue)
xlim([0 L+1])
xlabel('Residue')
ylabel('No. of amino acids')
post_proc_fig
mean_no_aa_at_a_site = mean(no_aa_site);
median_no_aa_at_a_site = median(no_aa_site);

%% Assigning weights to each seq w.r.t similar patient

seqs_noPatID = [];
for kk = 1:M    
    if strcmp(patientIDs{kk},'_')==1 %no patient ID info
        seqs_noPatID = [seqs_noPatID kk];
    end  
end

seqs_withPatID = setdiff(1:M,seqs_noPatID);

[unique_PatIDs, b,c] = unique(patientIDs(seqs_withPatID));

for kk = 1:length(c)
    w_patient(kk) = 1/length(find(c==c(kk)));  
end

%Assigning weight to each sequence
w = zeros(M,1);
w(seqs_noPatID) = 1; %Assigning a weight of 1 for seqs with no patient IDs
w(seqs_withPatID) = w_patient;

X = msa_aa;

N = sum(w);

%% Finding 100% conserved sites
[indices_conserved_sites,no_of_conserved_sites] = find_conserved_sites(X);

%% Statistics

fprintf('-----------------------------------------------------------------------------------\n')
fprintf('Sequence data (MSA) statistics\n')
fprintf('-----------------------------------------------------------------------------------\n')
fprintf('Number of residues, L = %d\n',L)
fprintf('Number of conserved residues = %d\n',no_of_conserved_sites)
fprintf('Number of mutating residues = %d\n', L - no_of_conserved_sites)
fprintf('Number of sequences = %d\n', M)
fprintf('Number of unique patients, N = %d\n',round(N))