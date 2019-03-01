function energy_vs_fitness(method, msa_aa_ex, msa_aa, phi_curr, mutant_order, ...
    H, ind_non_conserve)

% Code for comparing predicted energy with in vitro fitness measurements of
% E2 reported in literature
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%

run startup.m

% method = 1; %1-->including Js; 2-->only hs

% load test886_hepc_e2_send_Ahmed.mat
% 
% msa_aa_ex = msa_bin;
% phi_curr = num_mutants_combine_array;
% mutant_order = amino_single_combine_array;

if method == 1
%     H = reshape(J_mat,length(single_double_mutant_mat),length(single_double_mutant_mat));
    H = convertRayFLParamsToJohnParamsFormat(H);
%     H = triu(H);
elseif method == 2
    pi = mean(msa_aa_ex);
    hi_est_1pt = log(pi./(1-pi));
    H = diag(hi_est_1pt);
end

phi_cumulative(1) = phi_curr(1);
for kk = 2:length(mutant_order)
    phi_cumulative(kk) = sum(phi_curr(1:kk));
end

%Generating bin_matrix (definition of binary encoding)
for aba=1:30
    temp_matrix=[];
    temp_matrix = fliplr(eye(aba));
    temp_matrix = [zeros(1,size(temp_matrix,2)) ; temp_matrix ];
    bin_matrix{aba} = temp_matrix;
end

%% H77 sequence (E2 in H77 seq -- 384:746)

H77 = msa_aa(2,:);
H77_mut = H77(ind_non_conserve);

%% Paper 1:[Gal-Tanamy2008]


H77_mut_E655G = H77_mut;
H77_mut_E655G(rev_translation_indices(655-383,ind_non_conserve));
mutant_order{rev_translation_indices(655-383,ind_non_conserve)};
H77_mut_E655G(rev_translation_indices(655-383,ind_non_conserve)) = 'A';

H77_mut_N415Y = H77_mut;
H77_mut_N415Y(rev_translation_indices(415-383,ind_non_conserve));
mutant_order{rev_translation_indices(415-383,ind_non_conserve)};
H77_mut_N415Y(rev_translation_indices(415-383,ind_non_conserve)) = 'D';

H77_mut_E655G_N415Y = H77_mut;
H77_mut_E655G_N415Y(rev_translation_indices(415-383,ind_non_conserve)) = 'D';
H77_mut_E655G_N415Y(rev_translation_indices(655-383,ind_non_conserve)) = 'A';


msaToTest = [H77_mut;H77_mut_E655G;H77_mut_N415Y;H77_mut_E655G_N415Y];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq1 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq1(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

fitness1 = [800   900    70    20];
FFU1 = fitness1./fitness1(1);

FFU1_norm = standardize_data(fitness1);
Energy1_norm = standardize_data(EnergySeq1);

[r1,p1] = corr(Energy1_norm.',FFU1_norm.','type','spearman');


%% Paper 2:[Keck2009] H77 mutants

H77_mut_Q444A = H77_mut;
H77_mut_Q444A(rev_translation_indices(444-383,ind_non_conserve));
mutant_order{rev_translation_indices(444-383,ind_non_conserve)};
H77_mut_Q444A(rev_translation_indices(444-383,ind_non_conserve)) = 'A';

H77_mut_K446A = H77_mut;
H77_mut_K446A(rev_translation_indices(446-383,ind_non_conserve));
mutant_order{rev_translation_indices(446-383,ind_non_conserve)};
H77_mut_K446A(rev_translation_indices(446-383,ind_non_conserve)) = 'A';

H77_mut_E482A = H77_mut;
H77_mut_E482A(rev_translation_indices(482-383,ind_non_conserve));
mutant_order{rev_translation_indices(482-383,ind_non_conserve)};
H77_mut_E482A(rev_translation_indices(482-383,ind_non_conserve)) = 'Q';

H77_mut_S501A = H77_mut;
H77_mut_S501A(rev_translation_indices(501-383,ind_non_conserve));
mutant_order{rev_translation_indices(501-383,ind_non_conserve)};
H77_mut_S501A(rev_translation_indices(501-383,ind_non_conserve)) = 'A';

H77_mut_V506A = H77_mut;
H77_mut_V506A(rev_translation_indices(506-383,ind_non_conserve));
mutant_order{rev_translation_indices(506-383,ind_non_conserve)};
H77_mut_V506A(rev_translation_indices(506-383,ind_non_conserve)) = 'A';


msaToTest = [H77_mut;H77_mut_Q444A;H77_mut_K446A;H77_mut_E482A;H77_mut_S501A;H77_mut_V506A];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq2 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq2(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

fitness2 = [1700000     1000000      400000     2000000      700000       25000];
FFU2 = fitness2./fitness2(1);
FFU2_norm = standardize_data(fitness2);
Energy2_norm = standardize_data(EnergySeq2);

[r2,p2] = corr(Energy2_norm.',FFU2_norm.','type','spearman');



%% Paper 2b:[Keck2009] 02E10 mutants

[h,E10] = fastaread('H02.E10.fasta'); %synthetic sequence H02.E10
E10 = E10(end-363+1:end);
E10_mut = E10(ind_non_conserve);


E10_mut_Y444H = E10_mut;
E10_mut_Y444H(rev_translation_indices(444-383,ind_non_conserve));
mutant_order{rev_translation_indices(444-383,ind_non_conserve)};
E10_mut_Y444H(rev_translation_indices(444-383,ind_non_conserve)) = 'H';

E10_mut_Y444Q = E10_mut;
E10_mut_Y444Q(rev_translation_indices(444-383,ind_non_conserve));
mutant_order{rev_translation_indices(444-383,ind_non_conserve)};
E10_mut_Y444Q(rev_translation_indices(444-383,ind_non_conserve)) = 'V';

E10_mut_R446K = E10_mut;
E10_mut_R446K(rev_translation_indices(446-383,ind_non_conserve));
mutant_order{rev_translation_indices(446-383,ind_non_conserve)};
E10_mut_R446K(rev_translation_indices(446-383,ind_non_conserve)) = 'K';

E10_mut_R446G = E10_mut;
E10_mut_R446G(rev_translation_indices(446-383,ind_non_conserve));
mutant_order{rev_translation_indices(446-383,ind_non_conserve)};
E10_mut_R446G(rev_translation_indices(446-383,ind_non_conserve)) = 'Q';

E10_mut_Q482E = E10_mut;
E10_mut_Q482E(rev_translation_indices(482-383,ind_non_conserve));
mutant_order{rev_translation_indices(482-383,ind_non_conserve)};
E10_mut_Q482E(rev_translation_indices(482-383,ind_non_conserve)) = 'E';

E10_mut_N501S = E10_mut;
E10_mut_N501S(rev_translation_indices(501-383,ind_non_conserve));
mutant_order{rev_translation_indices(501-383,ind_non_conserve)};
E10_mut_N501S(rev_translation_indices(501-383,ind_non_conserve)) = 'S';

E10_mut_A506V = E10_mut;
E10_mut_A506V(rev_translation_indices(506-383,ind_non_conserve));
mutant_order{rev_translation_indices(506-383,ind_non_conserve)};
E10_mut_A506V(rev_translation_indices(506-383,ind_non_conserve)) = 'V';

E10_mut_N501S_A506V = E10_mut;
E10_mut_N501S_A506V(rev_translation_indices(501-383,ind_non_conserve)) = 'S';
E10_mut_N501S_A506V(rev_translation_indices(506-383,ind_non_conserve)) = 'V';

msaToTest = [H77_mut;E10_mut;E10_mut_Y444H;E10_mut_Y444Q;E10_mut_R446K;...
    E10_mut_R446G;E10_mut_Q482E;E10_mut_N501S;E10_mut_A506V;E10_mut_N501S_A506V];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq2b = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq2b(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

fitness2b = [1700000 16000 17000 15500 16500 20000 21000 21500 800000 1600000];
mean_small_2b = mean(fitness2b(2:8));
fitness2b_approx = [1700000 mean_small_2b 800000 1600000];
EnergySeq2b_approx = [EnergySeq2b(1) mean(EnergySeq2b(2:8)) EnergySeq2b(end-1) EnergySeq2b(end)];

% EnergySeq2b = EnergySeq2b(2:end);
FFU2b = fitness2b./fitness2b(1);
FFU2b_approx = fitness2b_approx./fitness2b_approx(1);

% [r2b,p2b] = corr([EnergySeq2b-EnergySeq2b(1)].',FFU2b.','type','spearman')

FFU2b_approx_norm = standardize_data(fitness2b_approx);
Energy2b_approx_norm = standardize_data(EnergySeq2b_approx);

[r2b_approx,p2b_approx] = corr(Energy2b_approx_norm.',FFU2b_approx_norm.','type','spearman');

%% Paper 3: [Pierce2016]

[~,~,dataInput] = xlsread('Infectivity_E2data_Pierce2016.xlsx');

for kk = 1:size(dataInput,1)-1
    strainName3{kk} = dataInput{kk+1,1};
    mut_site3(kk) = str2double(dataInput{kk+1,1}(end-2:end));
    amino_mut3(kk) =  strainName3{kk}(1);
    fitness3(kk) = dataInput{kk+1,2};
    significant_value(kk) = dataInput{kk+1,4};
end

for kk = 2:length(fitness3)
    H77_mut_Pierce2016(kk,:) = H77_mut;
    H77_mut_Pierce2016(kk,rev_translation_indices(mut_site3(kk)-383,ind_non_conserve)) = amino_mut3(kk);
end



msaToTest = [H77_mut;H77_mut_Pierce2016(2:end,:)];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq3 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq3(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

values_used = [1 find(fitness3>=3e3 & significant_value==115)];


fitness33 = [fitness3(values_used)];
EnergySeq3 = [EnergySeq3(values_used)];

FFU3 = fitness33./fitness33(1); %significant (*) and greater than detection threshold 3e3

EnergySeq3(8) = max(abs(diag(H))); %EnergySeq3(1)-min(diag(H)); %mutation A403 doesnt exist 


FFU3_norm = standardize_data(fitness33);
Energy3_norm = standardize_data(EnergySeq3);

[r3,p3] = corr(Energy3_norm.',FFU3_norm.','type','spearman');

%% Paper 4 [Goffard2005]

H77_mut_N417Q = H77_mut;
H77_mut_N417Q(rev_translation_indices(417-383,ind_non_conserve));
mutant_order{rev_translation_indices(417-383,ind_non_conserve)};
H77_mut_N417Q(rev_translation_indices(417-383,ind_non_conserve)) = 'D';

H77_mut_N423Q = H77_mut;
H77_mut_N423Q(rev_translation_indices(423-383,ind_non_conserve));
mutant_order{rev_translation_indices(423-383,ind_non_conserve)};
H77_mut_N423Q(rev_translation_indices(423-383,ind_non_conserve)) = 'S';

H77_mut_N430Q = H77_mut;
H77_mut_N430Q(rev_translation_indices(430-383,ind_non_conserve));
mutant_order{rev_translation_indices(430-383,ind_non_conserve)};
H77_mut_N430Q(rev_translation_indices(430-383,ind_non_conserve)) = 'E';

H77_mut_N448Q = H77_mut;
H77_mut_N448Q(rev_translation_indices(448-383,ind_non_conserve));
mutant_order{rev_translation_indices(448-383,ind_non_conserve)};
H77_mut_N448Q(rev_translation_indices(448-383,ind_non_conserve)) = 'D';

H77_mut_N476Q = H77_mut;
H77_mut_N476Q(rev_translation_indices(476-383,ind_non_conserve));
mutant_order{rev_translation_indices(476-383,ind_non_conserve)};
H77_mut_N476Q(rev_translation_indices(476-383,ind_non_conserve))= 'G';

H77_mut_N532Q = H77_mut;
H77_mut_N532Q(rev_translation_indices(532-383,ind_non_conserve));
mutant_order{rev_translation_indices(532-383,ind_non_conserve)};
H77_mut_N532Q(rev_translation_indices(532-383,ind_non_conserve)) = 'D';

H77_mut_N540Q = H77_mut;
H77_mut_N540Q(rev_translation_indices(540-383,ind_non_conserve));
mutant_order{rev_translation_indices(540-383,ind_non_conserve)};
H77_mut_N540Q(rev_translation_indices(540-383,ind_non_conserve)) = 'T';

H77_mut_N556Q = H77_mut;
H77_mut_N556Q(rev_translation_indices(556-383,ind_non_conserve));
mutant_order{rev_translation_indices(556-383,ind_non_conserve)};
H77_mut_N556Q(rev_translation_indices(556-383,ind_non_conserve)) = 'D';

H77_mut_N576Q = H77_mut;
H77_mut_N576Q(rev_translation_indices(576-383,ind_non_conserve));
mutant_order{rev_translation_indices(576-383,ind_non_conserve)};
H77_mut_N576Q(rev_translation_indices(576-383,ind_non_conserve)) = 'Q';

H77_mut_N623Q = H77_mut;
H77_mut_N623Q(rev_translation_indices(623-383,ind_non_conserve));
mutant_order{rev_translation_indices(623-383,ind_non_conserve)};
H77_mut_N623Q(rev_translation_indices(623-383,ind_non_conserve)) = 'Q';

H77_mut_N645Q = H77_mut;
H77_mut_N645Q(rev_translation_indices(645-383,ind_non_conserve));
mutant_order{rev_translation_indices(645-383,ind_non_conserve)};
H77_mut_N645Q(rev_translation_indices(645-383,ind_non_conserve)) = 'R';


msaToTest = [H77_mut;H77_mut_N417Q;H77_mut_N423Q;H77_mut_N430Q;H77_mut_N448Q;...
    H77_mut_N476Q;H77_mut_N532Q;H77_mut_N540Q;H77_mut_N556Q;H77_mut_N576Q;...
    H77_mut_N623Q;H77_mut_N645Q];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq6 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq6(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

fitness6 = [100	43.50282486	1.694915254	98.8700565	5.649717514	71.18644068	67.79661017	98.30508475	2.259887006	90.39548023	1.694915254	27.12];
FFU6 = fitness6./fitness6(1);

FFU6_norm = standardize_data(fitness6);
Energy6_norm = standardize_data(EnergySeq6);

[r6,p6] = corr(Energy6_norm.',FFU6_norm.','type','spearman');


%% Paper 7 [Falkowska2007]

H77_mut_T385A = H77_mut;
H77_mut_T385A(rev_translation_indices(385-383,ind_non_conserve));
mutant_order{rev_translation_indices(385-383,ind_non_conserve)};
H77_mut_T385A(rev_translation_indices(385-383,ind_non_conserve)) = 'A';

H77_mut_T388A = H77_mut;
H77_mut_T388A(rev_translation_indices(388-383,ind_non_conserve));
mutant_order{rev_translation_indices(388-383,ind_non_conserve)};
H77_mut_T388A(rev_translation_indices(388-383,ind_non_conserve)) = 'I';

H77_mut_N417A = H77_mut;
H77_mut_N417A(rev_translation_indices(417-383,ind_non_conserve));
mutant_order{rev_translation_indices(417-383,ind_non_conserve)};
H77_mut_N417A(rev_translation_indices(417-383,ind_non_conserve)) = 'A';

H77_mut_N423A = H77_mut;
H77_mut_N423A(rev_translation_indices(423-383,ind_non_conserve));
mutant_order{rev_translation_indices(423-383,ind_non_conserve)};
H77_mut_N423A(rev_translation_indices(423-383,ind_non_conserve)) = 'A';

H77_mut_N430A = H77_mut;
H77_mut_N430A(rev_translation_indices(430-383,ind_non_conserve));
mutant_order{rev_translation_indices(430-383,ind_non_conserve)};
H77_mut_N430A(rev_translation_indices(430-383,ind_non_conserve)) = 'A';

H77_mut_N448A = H77_mut;
H77_mut_N448A(rev_translation_indices(448-383,ind_non_conserve));
mutant_order{rev_translation_indices(448-383,ind_non_conserve)};
H77_mut_N448A(rev_translation_indices(448-383,ind_non_conserve)) = 'A';

H77_mut_N476A = H77_mut;
H77_mut_N476A(rev_translation_indices(476-383,ind_non_conserve));
mutant_order{rev_translation_indices(476-383,ind_non_conserve)};
H77_mut_N476A(rev_translation_indices(476-383,ind_non_conserve)) = 'A';

H77_mut_N510A = H77_mut;
H77_mut_N510A(rev_translation_indices(510-383,ind_non_conserve));
mutant_order{rev_translation_indices(510-383,ind_non_conserve)};
H77_mut_N510A(rev_translation_indices(510-383,ind_non_conserve)) = 'A';

H77_mut_N519A = H77_mut;
H77_mut_N519A(rev_translation_indices(519-383,ind_non_conserve));
mutant_order{rev_translation_indices(519-383,ind_non_conserve)};
H77_mut_N519A(rev_translation_indices(519-383,ind_non_conserve)) = 'A';

H77_mut_N532A = H77_mut;
H77_mut_N532A(rev_translation_indices(532-383,ind_non_conserve));
mutant_order{rev_translation_indices(532-383,ind_non_conserve)};
H77_mut_N532A(rev_translation_indices(532-383,ind_non_conserve)) = 'D';

H77_mut_N540A = H77_mut;
H77_mut_N540A(rev_translation_indices(540-383,ind_non_conserve));
mutant_order{rev_translation_indices(540-383,ind_non_conserve)};
H77_mut_N540A(rev_translation_indices(540-383,ind_non_conserve)) = 'A';

H77_mut_N556A = H77_mut;
H77_mut_N556A(rev_translation_indices(556-383,ind_non_conserve));
mutant_order{rev_translation_indices(556-383,ind_non_conserve)};
H77_mut_N556A(rev_translation_indices(556-383,ind_non_conserve)) = 'A';

H77_mut_N576A = H77_mut;
H77_mut_N576A(rev_translation_indices(576-383,ind_non_conserve));
mutant_order{rev_translation_indices(576-383,ind_non_conserve)};
H77_mut_N576A(rev_translation_indices(576-383,ind_non_conserve)) = 'A';

H77_mut_N623A = H77_mut;
H77_mut_N623A(rev_translation_indices(623-383,ind_non_conserve));
mutant_order{rev_translation_indices(623-383,ind_non_conserve)};
H77_mut_N623A(rev_translation_indices(623-383,ind_non_conserve)) = 'A';

H77_mut_N645A = H77_mut;
H77_mut_N645A(rev_translation_indices(645-383,ind_non_conserve));
mutant_order{rev_translation_indices(645-383,ind_non_conserve)};
H77_mut_N645A(rev_translation_indices(645-383,ind_non_conserve)) = 'A';


msaToTest = [H77_mut;H77_mut_T385A;H77_mut_T388A;H77_mut_N417A;H77_mut_N423A;...
    H77_mut_N430A;H77_mut_N448A;H77_mut_N476A;H77_mut_N510A;H77_mut_N519A;...
    H77_mut_N532A;H77_mut_N540A;H77_mut_N556A;H77_mut_N576A;H77_mut_N623A;H77_mut_N645A];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq7 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq7(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

indices_remove = [];%[4 6 11]+1; %Mutations not observed in data

fitness7 = [100 21.68674699	8.036144578	9.44578313	5.421686747	61.44578313	...
    6.024096386	94.57831325	2.409638554	3.012048193	33.73493976	85.54216867	...
    1.204819277	88.55421687	9.84337349 4.81927710843373];


FFU7 = fitness7./fitness7(1);

indices_small = find(fitness7<=10);
fitness7_approx = [fitness7(setdiff(1:length(fitness7),indices_small)) mean(fitness7(indices_small))];
EnergySeq7_approx = [EnergySeq7(setdiff(1:length(fitness7),indices_small)) mean(EnergySeq7(indices_small))];
FFU7_approx=fitness7_approx./fitness7_approx(1);

[r7_approx,p7_approx] = corr([EnergySeq7_approx-EnergySeq7_approx(1)].',FFU7_approx.','type','spearman');

%% Paper 8 [Rothwangl2008] - Huh7

H77_mut_Y474A = H77_mut;
H77_mut_Y474A(rev_translation_indices(474-383,ind_non_conserve));
mutant_order{rev_translation_indices(474-383,ind_non_conserve)};
H77_mut_Y474A(rev_translation_indices(474-383,ind_non_conserve)) = 'A';

H77_mut_D481A = H77_mut;
H77_mut_D481A(rev_translation_indices(481-383,ind_non_conserve));
mutant_order{rev_translation_indices(481-383,ind_non_conserve)};
H77_mut_D481A(rev_translation_indices(481-383,ind_non_conserve)) = 'A';

H77_mut_R483A = H77_mut;
H77_mut_R483A(rev_translation_indices(483-383,ind_non_conserve));
mutant_order{rev_translation_indices(483-383,ind_non_conserve)};
H77_mut_R483A(rev_translation_indices(483-383,ind_non_conserve)) = 'A';

H77_mut_Y485A = H77_mut;
H77_mut_Y485A(rev_translation_indices(485-383,ind_non_conserve));
mutant_order{rev_translation_indices(485-383,ind_non_conserve)};
H77_mut_Y485A(rev_translation_indices(485-383,ind_non_conserve)) = 'A';

H77_mut_W487A = H77_mut;
H77_mut_W487A(rev_translation_indices(487-383,ind_non_conserve));
mutant_order{rev_translation_indices(487-383,ind_non_conserve)};
H77_mut_W487A(rev_translation_indices(487-383,ind_non_conserve)) = 'A';

H77_mut_H488A = H77_mut;
H77_mut_H488A(rev_translation_indices(488-383,ind_non_conserve));
mutant_order{rev_translation_indices(488-383,ind_non_conserve)};
H77_mut_H488A(rev_translation_indices(488-383,ind_non_conserve)) = 'A';

H77_mut_Y489A = H77_mut;
H77_mut_Y489A(rev_translation_indices(489-383,ind_non_conserve));
mutant_order{rev_translation_indices(489-383,ind_non_conserve)};
H77_mut_Y489A(rev_translation_indices(489-383,ind_non_conserve)) = 'A';

H77_mut_R492A = H77_mut;
H77_mut_R492A(rev_translation_indices(492-383,ind_non_conserve));
mutant_order{rev_translation_indices(492-383,ind_non_conserve)};
H77_mut_R492A(rev_translation_indices(492-383,ind_non_conserve)) = 'A'; %'A' doesnt matter

H77_mut_Y527A = H77_mut;
H77_mut_Y527A(rev_translation_indices(527-383,ind_non_conserve));
mutant_order{rev_translation_indices(527-383,ind_non_conserve)};
H77_mut_Y527A(rev_translation_indices(527-383,ind_non_conserve)) = 'A';

H77_mut_W529A = H77_mut;
H77_mut_W529A(rev_translation_indices(529-383,ind_non_conserve));
mutant_order{rev_translation_indices(529-383,ind_non_conserve)};
H77_mut_W529A(rev_translation_indices(529-383,ind_non_conserve)) = 'A';

H77_mut_D533A = H77_mut;
H77_mut_D533A(rev_translation_indices(533-383,ind_non_conserve));
mutant_order{rev_translation_indices(533-383,ind_non_conserve)};
H77_mut_D533A(rev_translation_indices(533-383,ind_non_conserve)) = 'A';

H77_mut_D535A = H77_mut;
H77_mut_D535A(rev_translation_indices(535-383,ind_non_conserve));
mutant_order{rev_translation_indices(535-383,ind_non_conserve)};
H77_mut_D535A(rev_translation_indices(535-383,ind_non_conserve)) = 'A';

H77_mut_W549A = H77_mut;
H77_mut_W549A(rev_translation_indices(549-383,ind_non_conserve));
mutant_order{rev_translation_indices(549-383,ind_non_conserve)};
H77_mut_W549A(rev_translation_indices(549-383,ind_non_conserve)) = 'A';

H77_mut_F550A = H77_mut;
H77_mut_F550A(rev_translation_indices(550-383,ind_non_conserve));
mutant_order{rev_translation_indices(550-383,ind_non_conserve)};
H77_mut_F550A(rev_translation_indices(550-383,ind_non_conserve)) = 'A';

H77_mut_Y613A = H77_mut;
H77_mut_Y613A(rev_translation_indices(613-383,ind_non_conserve));
mutant_order{rev_translation_indices(613-383,ind_non_conserve)};
H77_mut_Y613A(rev_translation_indices(613-383,ind_non_conserve)) = 'A';

H77_mut_R614A = H77_mut;
H77_mut_R614A(rev_translation_indices(614-383,ind_non_conserve));
mutant_order{rev_translation_indices(614-383,ind_non_conserve)};
H77_mut_R614A(rev_translation_indices(614-383,ind_non_conserve)) = 'A';

H77_mut_L615A = H77_mut;
H77_mut_L615A(rev_translation_indices(615-383,ind_non_conserve));
mutant_order{rev_translation_indices(615-383,ind_non_conserve)};
H77_mut_L615A(rev_translation_indices(615-383,ind_non_conserve)) = 'A';

H77_mut_W616A = H77_mut;
H77_mut_W616A(rev_translation_indices(616-383,ind_non_conserve));
mutant_order{rev_translation_indices(616-383,ind_non_conserve)};
H77_mut_W616A(rev_translation_indices(616-383,ind_non_conserve)) = 'A';

H77_mut_H617A = H77_mut;
H77_mut_H617A(rev_translation_indices(617-383,ind_non_conserve));
mutant_order{rev_translation_indices(617-383,ind_non_conserve)};
H77_mut_H617A(rev_translation_indices(617-383,ind_non_conserve)) = 'A';

H77_mut_Y618A = H77_mut;
H77_mut_Y618A(rev_translation_indices(618-383,ind_non_conserve));
mutant_order{rev_translation_indices(618-383,ind_non_conserve)};
H77_mut_Y618A(rev_translation_indices(618-383,ind_non_conserve)) = 'A';


msaToTest = [H77_mut;H77_mut_Y474A;H77_mut_D481A;H77_mut_R483A;H77_mut_Y485A;...
    H77_mut_W487A;H77_mut_H488A;H77_mut_Y489A;H77_mut_R492A;H77_mut_Y527A;...
    H77_mut_W529A;H77_mut_D533A;H77_mut_D535A;H77_mut_W549A;H77_mut_F550A;H77_mut_Y613A;...
    H77_mut_R614A;H77_mut_L615A;H77_mut_W616A;H77_mut_H617A;H77_mut_Y618A];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq8 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq8(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end
% 
fitness8 = [100 77 64 2 2 1 2 4 5 5 1 12 2 2 19 3 3 44 2 2 2]; %Table 1 Huh7
fitness8_approx = [fitness8(1:3)  mean(fitness8([4:11 13 14 16 17 19:21])) fitness8([12 15 18])];
EnergySeq8_approx = [EnergySeq8(1:3)  mean(EnergySeq8([4:11 13 14 16 17 19:21])) EnergySeq8([12 15 18])];

FFU8 = fitness8./fitness8(1);
FFU8_approx = fitness8_approx./fitness8_approx(1);

FFU8_approx_norm = standardize_data(fitness8_approx);
Energy8_approx_norm = standardize_data(EnergySeq8_approx);

[r8_approx,p8_approx] = corr(Energy8_approx_norm.',FFU8_approx_norm.','type','spearman');

%% Paper 8b [Rothwangl2008] - Hep3b

H77_mut_Y474A = H77_mut;
H77_mut_Y474A(rev_translation_indices(474-383,ind_non_conserve));
mutant_order{rev_translation_indices(474-383,ind_non_conserve)};
H77_mut_Y474A(rev_translation_indices(474-383,ind_non_conserve)) = 'A';

H77_mut_D481A = H77_mut;
H77_mut_D481A(rev_translation_indices(481-383,ind_non_conserve));
mutant_order{rev_translation_indices(481-383,ind_non_conserve)};
H77_mut_D481A(rev_translation_indices(481-383,ind_non_conserve)) = 'A';

H77_mut_R483A = H77_mut;
H77_mut_R483A(rev_translation_indices(483-383,ind_non_conserve));
mutant_order{rev_translation_indices(483-383,ind_non_conserve)};
H77_mut_R483A(rev_translation_indices(483-383,ind_non_conserve)) = 'A';

H77_mut_Y485A = H77_mut;
H77_mut_Y485A(rev_translation_indices(485-383,ind_non_conserve));
mutant_order{rev_translation_indices(485-383,ind_non_conserve)};
H77_mut_Y485A(rev_translation_indices(485-383,ind_non_conserve)) = 'A';

H77_mut_W487A = H77_mut;
H77_mut_W487A(rev_translation_indices(487-383,ind_non_conserve));
mutant_order{rev_translation_indices(487-383,ind_non_conserve)};
H77_mut_W487A(rev_translation_indices(487-383,ind_non_conserve)) = 'A';

H77_mut_H488A = H77_mut;
H77_mut_H488A(rev_translation_indices(488-383,ind_non_conserve));
mutant_order{rev_translation_indices(488-383,ind_non_conserve)};
H77_mut_H488A(rev_translation_indices(488-383,ind_non_conserve)) = 'A';

H77_mut_Y489A = H77_mut;
H77_mut_Y489A(rev_translation_indices(489-383,ind_non_conserve));
mutant_order{rev_translation_indices(489-383,ind_non_conserve)};
H77_mut_Y489A(rev_translation_indices(489-383,ind_non_conserve)) = 'A';

H77_mut_R492A = H77_mut;
H77_mut_R492A(rev_translation_indices(492-383,ind_non_conserve));
mutant_order{rev_translation_indices(492-383,ind_non_conserve)};
H77_mut_R492A(rev_translation_indices(492-383,ind_non_conserve)) = 'A';

H77_mut_Y527A = H77_mut;
H77_mut_Y527A(rev_translation_indices(527-383,ind_non_conserve));
mutant_order{rev_translation_indices(527-383,ind_non_conserve)};
H77_mut_Y527A(rev_translation_indices(527-383,ind_non_conserve)) = 'A';

H77_mut_W529A = H77_mut;
H77_mut_W529A(rev_translation_indices(529-383,ind_non_conserve));
mutant_order{rev_translation_indices(529-383,ind_non_conserve)};
H77_mut_W529A(rev_translation_indices(529-383,ind_non_conserve)) = 'A';

H77_mut_D533A = H77_mut;
H77_mut_D533A(rev_translation_indices(533-383,ind_non_conserve));
mutant_order{rev_translation_indices(533-383,ind_non_conserve)};
H77_mut_D533A(rev_translation_indices(533-383,ind_non_conserve)) = 'A';

H77_mut_D535A = H77_mut;
H77_mut_D535A(rev_translation_indices(535-383,ind_non_conserve));
mutant_order{rev_translation_indices(535-383,ind_non_conserve)};
H77_mut_D535A(rev_translation_indices(535-383,ind_non_conserve)) = 'A';

H77_mut_W549A = H77_mut;
H77_mut_W549A(rev_translation_indices(549-383,ind_non_conserve));
mutant_order{rev_translation_indices(549-383,ind_non_conserve)};
H77_mut_W549A(rev_translation_indices(549-383,ind_non_conserve)) = 'A';

H77_mut_F550A = H77_mut;
H77_mut_F550A(rev_translation_indices(550-383,ind_non_conserve));
mutant_order{rev_translation_indices(550-383,ind_non_conserve)};
H77_mut_F550A(rev_translation_indices(550-383,ind_non_conserve)) = 'A';

H77_mut_Y613A = H77_mut;
H77_mut_Y613A(rev_translation_indices(613-383,ind_non_conserve));
mutant_order{rev_translation_indices(613-383,ind_non_conserve)};
H77_mut_Y613A(rev_translation_indices(613-383,ind_non_conserve)) = 'A';

H77_mut_R614A = H77_mut;
H77_mut_R614A(rev_translation_indices(614-383,ind_non_conserve));
mutant_order{rev_translation_indices(614-383,ind_non_conserve)};
H77_mut_R614A(rev_translation_indices(614-383,ind_non_conserve)) = 'A';

H77_mut_L615A = H77_mut;
H77_mut_L615A(rev_translation_indices(615-383,ind_non_conserve));
mutant_order{rev_translation_indices(615-383,ind_non_conserve)};
H77_mut_L615A(rev_translation_indices(615-383,ind_non_conserve)) = 'A';

H77_mut_W616A = H77_mut;
H77_mut_W616A(rev_translation_indices(616-383,ind_non_conserve));
mutant_order{rev_translation_indices(616-383,ind_non_conserve)};
H77_mut_W616A(rev_translation_indices(616-383,ind_non_conserve)) = 'A';

H77_mut_H617A = H77_mut;
H77_mut_H617A(rev_translation_indices(617-383,ind_non_conserve));
mutant_order{rev_translation_indices(617-383,ind_non_conserve)};
H77_mut_H617A(rev_translation_indices(617-383,ind_non_conserve)) = 'A';

H77_mut_Y618A = H77_mut;
H77_mut_Y618A(rev_translation_indices(618-383,ind_non_conserve));
mutant_order{rev_translation_indices(618-383,ind_non_conserve)};
H77_mut_Y618A(rev_translation_indices(618-383,ind_non_conserve)) = 'A';


msaToTest = [H77_mut;H77_mut_Y474A;H77_mut_D481A;H77_mut_R483A;H77_mut_Y485A;...
    H77_mut_W487A;H77_mut_H488A;H77_mut_Y489A;H77_mut_R492A;H77_mut_Y527A;...
    H77_mut_W529A;H77_mut_D533A;H77_mut_D535A;H77_mut_W549A;H77_mut_F550A;H77_mut_Y613A;...
    H77_mut_R614A;H77_mut_L615A;H77_mut_W616A;H77_mut_H617A;H77_mut_Y618A];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq8b = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq8b(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

fitness8b = [100 36 71 1 1 0 1 2 2 1 1 5 0 0 2 1 0 10 0 0 1]; %Table 1 Hep3B

FFU8b = fitness8b./fitness8b(1);

indices_small = find(fitness8b<10);
fitness8b_approx = [fitness8b(setdiff(1:length(fitness8b),indices_small)) mean(fitness8b(indices_small))];
EnergySeq8b_approx = [EnergySeq8b(setdiff(1:length(fitness8b),indices_small)) mean(EnergySeq8b(indices_small))];

FFU8b_approx = fitness8b_approx./fitness8b_approx(1);

[r8b_approx, p8b_approx] = corr([EnergySeq8b_approx-EnergySeq8b_approx(1)].',FFU8b_approx.','type','spearman');


%% Paper 9: [Drummer2006]

H77_mut_G436A = H77_mut;
H77_mut_G436A(rev_translation_indices(436-383,ind_non_conserve));
mutant_order{rev_translation_indices(436-383,ind_non_conserve)};
H77_mut_G436A(rev_translation_indices(436-383,ind_non_conserve)) = 'D';

H77_mut_W437F = H77_mut;
H77_mut_W437F(rev_translation_indices(437-383,ind_non_conserve));
mutant_order{rev_translation_indices(437-383,ind_non_conserve)};
H77_mut_W437F(rev_translation_indices(437-383,ind_non_conserve)) = 'F';

H77_mut_L438V = H77_mut;
H77_mut_L438V(rev_translation_indices(438-383,ind_non_conserve));
mutant_order{rev_translation_indices(438-383,ind_non_conserve)};
H77_mut_L438V(rev_translation_indices(438-383,ind_non_conserve)) = 'V';

H77_mut_A439G = H77_mut; %*
H77_mut_A439G(rev_translation_indices(439-383,ind_non_conserve));
mutant_order{rev_translation_indices(439-383,ind_non_conserve)};
H77_mut_A439G(rev_translation_indices(439-383,ind_non_conserve)) = 'T';

H77_mut_G440A = H77_mut;
H77_mut_G440A(rev_translation_indices(440-383,ind_non_conserve));
mutant_order{rev_translation_indices(440-383,ind_non_conserve)};
H77_mut_G440A(rev_translation_indices(440-383,ind_non_conserve)) = 'A';

H77_mut_L441V = H77_mut; %*
H77_mut_L441V(rev_translation_indices(441-383,ind_non_conserve));
mutant_order{rev_translation_indices(441-383,ind_non_conserve)};
H77_mut_L441V(rev_translation_indices(441-383,ind_non_conserve)) = 'P'; 

H77_mut_F442L = H77_mut;
H77_mut_F442L(rev_translation_indices(442-383,ind_non_conserve));
mutant_order{rev_translation_indices(442-383,ind_non_conserve)};
H77_mut_F442L(rev_translation_indices(442-383,ind_non_conserve)) = 'L';

H77_mut_F442M = H77_mut;
H77_mut_F442M(rev_translation_indices(442-383,ind_non_conserve));
mutant_order{rev_translation_indices(442-383,ind_non_conserve)};
H77_mut_F442M(rev_translation_indices(442-383,ind_non_conserve)) = 'I'; 

H77_mut_Y443A = H77_mut; %*
H77_mut_Y443A(rev_translation_indices(443-383,ind_non_conserve));
mutant_order{rev_translation_indices(443-383,ind_non_conserve)};
H77_mut_Y443A(rev_translation_indices(443-383,ind_non_conserve)) = 'S';

msaToTest = [H77_mut;H77_mut_G436A;H77_mut_W437F;H77_mut_L438V;H77_mut_A439G;...
    H77_mut_G440A;H77_mut_L441V;H77_mut_F442L;H77_mut_F442M;H77_mut_Y443A];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq9 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq9(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

fitness9 = [100 10.31746032	40.21164021	12.6984127	9.469 19.57671958	4.637 62.96296296 42 6.878306878];

FFU9 = fitness9./fitness9(1);

[r9,p9] = corr([EnergySeq9-EnergySeq9(1)].',FFU9.','type','spearman');


%% Paper10 [Keck2012] Fig. 6A

H77_mut_W420A = H77_mut;
H77_mut_W420A(rev_translation_indices(420-383,ind_non_conserve));
mutant_order{rev_translation_indices(420-383,ind_non_conserve)};
H77_mut_W420A(rev_translation_indices(420-383,ind_non_conserve)) = 'A';%R

H77_mut_N428A = H77_mut;
H77_mut_N428A(rev_translation_indices(428-383,ind_non_conserve));
mutant_order{rev_translation_indices(428-383,ind_non_conserve)};
H77_mut_N428A(rev_translation_indices(428-383,ind_non_conserve)) = 'A';%ST

H77_mut_W437A = H77_mut;
H77_mut_W437A(rev_translation_indices(437-383,ind_non_conserve));
mutant_order{rev_translation_indices(437-383,ind_non_conserve)};
H77_mut_W437A(rev_translation_indices(437-383,ind_non_conserve)) = 'A'; %F

H77_mut_L441A = H77_mut;
H77_mut_L441A(rev_translation_indices(441-383,ind_non_conserve));
mutant_order{rev_translation_indices(441-383,ind_non_conserve)};
H77_mut_L441A(rev_translation_indices(441-383,ind_non_conserve)) = 'A'; %P

H77_mut_F442A = H77_mut;
H77_mut_F442A(rev_translation_indices(442-383,ind_non_conserve));
mutant_order{rev_translation_indices(442-383,ind_non_conserve)};
H77_mut_F442A(rev_translation_indices(442-383,ind_non_conserve)) = 'A'; %LI

H77_mut_Y443A = H77_mut;
H77_mut_Y443A(rev_translation_indices(443-383,ind_non_conserve));
mutant_order{rev_translation_indices(443-383,ind_non_conserve)};
H77_mut_Y443A(rev_translation_indices(443-383,ind_non_conserve)) = 'A'; %HS

H77_mut_K446A = H77_mut;
H77_mut_K446A(rev_translation_indices(446-383,ind_non_conserve));
mutant_order{rev_translation_indices(446-383,ind_non_conserve)};
H77_mut_K446A(rev_translation_indices(446-383,ind_non_conserve)) = 'A'; %RNQ

H77_mut_Y613A = H77_mut;
H77_mut_Y613A(rev_translation_indices(613-383,ind_non_conserve));
mutant_order{rev_translation_indices(613-383,ind_non_conserve)};
H77_mut_Y613A(rev_translation_indices(613-383,ind_non_conserve)) = 'A'; %H

H77_mut_W616A = H77_mut;
H77_mut_W616A(rev_translation_indices(616-383,ind_non_conserve));
mutant_order{rev_translation_indices(616-383,ind_non_conserve)};
H77_mut_W616A(rev_translation_indices(616-383,ind_non_conserve)) = 'A'; %R


msaToTest = [H77_mut;H77_mut_W420A;H77_mut_N428A;H77_mut_W437A;H77_mut_L441A;H77_mut_F442A;...
    H77_mut_Y443A;H77_mut_K446A;H77_mut_Y613A;H77_mut_W616A];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq10 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq10(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

fitness10 = [34.074	2.751	7.619	2.751	2.751	2.751	5.926	13.333	2.328	2.328];
fitness10_approx = [fitness10([1 3 7 8]) mean(fitness10([2 4 5 6 9 10]))];
EnergySeq10_approx = [EnergySeq10([1 3 7 8]) mean(EnergySeq10([2 4 5 6 9 10]))];


FFU10 = fitness10./fitness10(1);
FFU10_approx = fitness10_approx./fitness10_approx(1);


FFU10_norm = standardize_data(fitness10);
Energy10_norm = standardize_data(EnergySeq10);

FFU10_approx_norm = standardize_data(fitness10_approx);
Energy10_approx_norm = standardize_data(EnergySeq10_approx);

[r10_approx,p10_approx] = corr(Energy10_approx_norm.',FFU10_approx_norm.','type','spearman');

% fitness_approx(FFU10,EnergySeq10)

%% Paper 11 [Gopal2017] Table S5

mutations = 600:645;
clear msaToTest
for kk = 1:length(mutations)
    mut_seq = H77_mut;
    mut_seq(rev_translation_indices(mutations(kk)-383,ind_non_conserve)) = 'A';
    msaToTest(kk,:) = mut_seq;
end

msaToTest = [H77_mut;msaToTest];

numSeqMsaToTest = size(msaToTest,1);
EnergySeq11 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq11(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end


fitness11 = [100 5	13	2	1	5	12	2	3	2	6	1	1	4	3	2	4 ...
    0	1	0	1	0	1	11	1	3	2	7	9	5	27	4	4	1	0	1 ...
    1	1	5	0	1	0	25	3	1	1	6];
FFU11 = fitness11./fitness11(1);

indices_exclude = find(fitness11<10);

fitness_small = mean(fitness11(indices_exclude));
EnergySeq11_small = mean(EnergySeq11(indices_exclude));

fitness11(indices_exclude) = [];
fitness11 = [fitness11 fitness_small];

EnergySeq11_approx = EnergySeq11;
EnergySeq11_approx(indices_exclude) = [];
EnergySeq11_approx = [EnergySeq11_approx EnergySeq11_small];

FFU11_approx = fitness11./fitness11(1);

[r11_approx,p11_approx] = corr([EnergySeq11_approx - EnergySeq11_approx(1)].',FFU11_approx.','type','spearman');


%% Paper 14 - El-Diwany2017 
%used S5FigA 
%   - more significant result ***
%   - 1a154 strain is H77

H77_mut_L403F = H77_mut;
H77_mut_L403F(rev_translation_indices(403-383,ind_non_conserve));
mutant_order{rev_translation_indices(403-383,ind_non_conserve)};
H77_mut_L403F(rev_translation_indices(403-383,ind_non_conserve)) = 'F';%R

H77_mut_L438V = H77_mut;
H77_mut_L438V(rev_translation_indices(438-383,ind_non_conserve));
mutant_order{rev_translation_indices(438-383,ind_non_conserve)};
H77_mut_L438V(rev_translation_indices(438-383,ind_non_conserve)) = 'V';%R

msaToTest = [H77_mut;H77_mut_L403F;H77_mut_L438V];

numSeqMsaToTest = size(msaToTest,1);
EnergySeq14 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq14(i)] = calcSeqEnergy(out_seq_ex(i,:),H);    
    
end

fitness14 = [1e7 8e6 1e6];
FFU14 = fitness14./fitness14(1);
[r14,p14] = corr([EnergySeq14-EnergySeq14(1)].',FFU14.','type','spearman');


%% Paper 15 - Guan2012

H77_mut_A397G = H77_mut;
H77_mut_A397G(rev_translation_indices(397-383,ind_non_conserve));
mutant_order{rev_translation_indices(397-383,ind_non_conserve)};
H77_mut_A397G(rev_translation_indices(397-383,ind_non_conserve)) = 'G';

H77_mut_G398A = H77_mut;
H77_mut_G398A(rev_translation_indices(398-383,ind_non_conserve));
mutant_order{rev_translation_indices(398-383,ind_non_conserve)};
H77_mut_G398A(rev_translation_indices(398-383,ind_non_conserve)) = 'A';

H77_mut_A397G_G398A = H77_mut;
H77_mut_A397G_G398A(rev_translation_indices([397 398]-383,ind_non_conserve));
mutant_order{rev_translation_indices(398-383,ind_non_conserve)};
H77_mut_A397G_G398A(rev_translation_indices([397 398]-383,ind_non_conserve)) = 'GA';

H77_mut_K408A = H77_mut;
H77_mut_K408A(rev_translation_indices(408-383,ind_non_conserve));
mutant_order{rev_translation_indices(408-383,ind_non_conserve)};
H77_mut_K408A(rev_translation_indices(408-383,ind_non_conserve)) = 'A';

H77_mut_Q409A = H77_mut;
H77_mut_Q409A(rev_translation_indices(409-383,ind_non_conserve));
mutant_order{rev_translation_indices(409-383,ind_non_conserve)};
H77_mut_Q409A(rev_translation_indices(409-383,ind_non_conserve)) = 'A';

H77_mut_N410A = H77_mut;
H77_mut_N410A(rev_translation_indices(410-383,ind_non_conserve));
mutant_order{rev_translation_indices(410-383,ind_non_conserve)};
H77_mut_N410A(rev_translation_indices(410-383,ind_non_conserve)) = 'A';

msaToTest = [H77_mut;H77_mut_A397G;H77_mut_G398A;H77_mut_A397G_G398A;...
    H77_mut_K408A;H77_mut_Q409A;H77_mut_N410A];

numSeqMsaToTest = size(msaToTest,1);
EnergySeq15 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq15(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
    
    
end

fitness15 = [99.31	79.31	73.2  39.31	80.01	73.103	72.793];

FFU15 = fitness15./fitness15(1);
[r15,p15] = corr([EnergySeq15-EnergySeq15(1)].',FFU15.','type','spearman');



%% Fitness vs Energy (normalized)

% set(0,'DefaultAxesFontName','Arial')
% set(0,'DefaultTextFontName','Arial')
% set(0,'DefaultAxesFontSize',6)
% set(0,'DefaultTextFontSize',6)

markersize = 6;
line_width = 0.5;

Energy1_norm = standardize_data(EnergySeq1);
Energy2_norm = standardize_data(EnergySeq2);
Energy2b_approx_norm = standardize_data(EnergySeq2b_approx);
Energy3_norm = standardize_data(EnergySeq3);
Energy6_norm = standardize_data(EnergySeq6);
Energy7_approx_norm = standardize_data(EnergySeq7_approx);
Energy8_approx_norm = standardize_data(EnergySeq8_approx);
Energy8b_approx_norm = standardize_data(EnergySeq8b_approx);
Energy9_norm = standardize_data(EnergySeq9);
Energy10_approx_norm = standardize_data(EnergySeq10_approx);
Energy11_approx_norm = standardize_data(EnergySeq11_approx);
Energy14_norm = standardize_data(EnergySeq14);
Energy15_norm = standardize_data(EnergySeq15);


FFU1_norm = standardize_data(FFU1);
FFU2_norm = standardize_data(FFU2);
FFU2b_approx_norm = standardize_data(FFU2b_approx);
FFU3_norm = standardize_data(FFU3);
FFU6_norm = standardize_data(FFU6);
FFU7_approx_norm = standardize_data(FFU7_approx);
FFU8_approx_norm = standardize_data(FFU8_approx);
FFU8b_approx_norm = standardize_data(FFU8b_approx);
FFU9_norm = standardize_data(FFU9);
FFU10_approx_norm = standardize_data(FFU10_approx);
FFU11_approx_norm = standardize_data(FFU11_approx);
FFU14_norm = standardize_data(FFU14);
FFU15_norm = standardize_data(FFU15);

pink = color_scheme(8,:);
gray = [0.2 0.2 0.2];
figure;
plot(Energy6_norm,FFU6_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',red,'LineWidth',line_width)

xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')

hold on 
plot(Energy9_norm,FFU9_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',blue,'LineWidth',line_width)
plot(Energy7_approx_norm,FFU7_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',green,'LineWidth',line_width)
plot(Energy8_approx_norm,FFU8_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'LineWidth',line_width)
plot(Energy8b_approx_norm,FFU8b_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','k','LineWidth',line_width)
plot(Energy1_norm,FFU1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',orange,'Linewidth',line_width)
% plot(Energy5,titer5,'<','MarkerEdgeColor',purple,'MarkerSize',10,'Linewidth',1.5)
plot(Energy2_norm,FFU2_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',yellow,'LineWidth',line_width)
plot(Energy2b_approx_norm,FFU2b_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',lightblue,'Linewidth',line_width)%,'MarkerFaceColor','r')
plot(Energy10_approx_norm,FFU10_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',pink,'Linewidth',line_width)
plot(Energy15_norm,FFU15_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','c','Linewidth',line_width)
plot(Energy3_norm,FFU3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',brown,'Linewidth',line_width)
plot(Energy11_approx_norm,FFU11_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',gray,'Linewidth',line_width)
plot(Energy14_norm,FFU14_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',lightgreen,'Linewidth',line_width)
% plot(Energy2,titer2,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',lightgreen,'Linewidth',line_width)
% legend('Gal-Tanamy2008','Keck2009-H77','Keck2009-02.E10','Pierce2016')
axis([-4 4 -3 5])

EnergyVecNorm = [Energy1_norm Energy2_norm Energy2b_approx_norm Energy3_norm Energy6_norm Energy7_approx_norm ...
    Energy8_approx_norm Energy8b_approx_norm Energy9_norm Energy10_approx_norm Energy11_approx_norm Energy14_norm];
FFUVecNorm = [FFU1_norm FFU2_norm FFU2b_approx_norm FFU3_norm FFU6_norm FFU7_approx_norm ...
    FFU8_approx_norm FFU8b_approx_norm FFU9_norm FFU10_approx_norm FFU11_approx_norm FFU14_norm];

P = polyfit(EnergyVecNorm',FFUVecNorm',1);
x = -4:5; %xaxis
y = P(1)*x+P(2);
plot(x,y,'k--','LineWidth',1)

FFU_vec = [FFU1 FFU2 FFU2b_approx FFU3 FFU6 FFU7_approx FFU8_approx FFU8b_approx FFU9 FFU10_approx FFU11_approx FFU14 FFU15];

total_length = length(FFU_vec);
w1 = length(FFU1)/total_length;
w2 = length(FFU2)/total_length;
w2b_approx = length(FFU2b_approx)/total_length;
w3 = length(FFU3)/total_length;
w6 = length(FFU6)/total_length;
w7_approx = length(FFU7_approx)/total_length;
w8_approx = length(FFU8_approx)/total_length;
w8b_approx = length(FFU8b_approx)/total_length;
w9 = length(FFU9)/total_length;
w10_approx = length(FFU10_approx)/total_length;
w11_approx = length(FFU11_approx)/total_length;
w14 = length(FFU14)/total_length;
w15 = length(FFU15)/total_length;


rho_weighted_average = sum([r1 r2 r2b_approx r3 r6 r7_approx r8_approx r8b_approx r9 r10_approx r11_approx r14 r15]...
    .*[w1 w2 w2b_approx w3 w6 w7_approx w8_approx w8b_approx w9 w10_approx w11_approx w14 w15]);

text(-3,-2,sprintf('$$\\bar{r} = %.2f$$',rho_weighted_average),'interpreter','latex','FontSize',16)

h = legend('(30)','(31)','(32)','(33)-Huh7','(33)-Hep3b',...
    '(34)','(16)','(16)-02.E10','(35)','(36)','(37)','(38)','(39)',...
    'LS fit','Location','NorthEast');

h.Position(1) = h.Position(1)+0.05;
h.Position(2) = h.Position(2)-0.1;

legend boxoff


set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.01 .01] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , -5:1:5, ...
  'XTick'       , -5:1:5, ...
  'LineWidth'   , 0.5        );





