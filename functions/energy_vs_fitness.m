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

% [~,~,dataInput] = xlsread('FFU_E2data_Gal-Tanamy2008.xlsx');
% 
% for kk = 1:size(dataInput,1)-1
%     strainName1{kk} = dataInput{kk+1,1};
%     fitness1(kk) = dataInput{kk+1,2};
% end

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

% figure;
% plot(Energy1_norm,FFU1_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r1,p1] = corr(Energy1_norm.',FFU1_norm.','type','spearman');
% % title(sprintf('Spearman r = %.3f, p = %.3e',r1,p1),'interpreter','tex')

%% Paper 2:[Keck2009] H77 mutants

% [~,~,dataInput] = xlsread('RLU_E2data_H77mutants_Keck2009.xlsx');
% 
% for kk = 1:size(dataInput,1)-1
%     strainName2{kk} = dataInput{kk+1,1};
%     fitness2(kk) = dataInput{kk+1,2};
% end

H77_mut_Q444A = H77_mut;
H77_mut_Q444A(rev_translation_indices(444-383,ind_non_conserve));
mutant_order{rev_translation_indices(444-383,ind_non_conserve)};
H77_mut_Q444A(rev_translation_indices(444-383,ind_non_conserve)) = 'A';

H77_mut_N415Y = H77_mut;
H77_mut_N415Y(rev_translation_indices(446-383,ind_non_conserve));
mutant_order{rev_translation_indices(446-383,ind_non_conserve)};
H77_mut_N415Y(rev_translation_indices(446-383,ind_non_conserve)) = 'Q';

H77_mut_E482A = H77_mut;
H77_mut_E482A(rev_translation_indices(482-383,ind_non_conserve));
mutant_order{rev_translation_indices(482-383,ind_non_conserve)};
H77_mut_E482A(rev_translation_indices(482-383,ind_non_conserve)) = 'Q';

H77_mut_S501A = H77_mut;
H77_mut_S501A(rev_translation_indices(501-383,ind_non_conserve));
mutant_order{rev_translation_indices(501-383,ind_non_conserve)};
H77_mut_S501A(rev_translation_indices(501-383,ind_non_conserve)) = 'K';

H77_mut_V506A = H77_mut;
H77_mut_V506A(rev_translation_indices(506-383,ind_non_conserve));
mutant_order{rev_translation_indices(506-383,ind_non_conserve)};
H77_mut_V506A(rev_translation_indices(506-383,ind_non_conserve)) = 'I';


msaToTest = [H77_mut;H77_mut_Q444A;H77_mut_N415Y;H77_mut_E482A;H77_mut_S501A;H77_mut_V506A];
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

% figure;
% plot([EnergySeq2-EnergySeq2(1)],FFU2,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r2,p2] = corr([EnergySeq2-EnergySeq2(1)].',FFU2.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r2,p2),'interpreter','tex')

FFU2_norm = standardize_data(fitness2);
Energy2_norm = standardize_data(EnergySeq2);

% figure;
% plot(Energy2_norm,FFU2_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r2,p2] = corr(Energy2_norm.',FFU2_norm.','type','spearman');
% title(sprintf('Spearman r = %.3f, p = %.3e',r2,p2),'interpreter','tex')

%% Paper 2b:[Keck2009] 02E10 mutants

[h,E10] = fastaread('H02.E10.fasta'); %synthetic sequence H02.E10
E10 = E10(end-363+1:end);
E10_mut = E10(ind_non_conserve);

% [~,~,dataInput] = xlsread('RLU_E2data_02E10mutants_Keck2009.xlsx');
% 
% for kk = 1:size(dataInput,1)-1
%     strainName2b{kk} = dataInput{kk+1,1};
%     fitness2b(kk) = dataInput{kk+1,2};
% end

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

msaToTest = [H77_mut;E10_mut;E10_mut_Y444H;E10_mut_Y444Q;E10_mut_R446K;E10_mut_R446G;E10_mut_Q482E;E10_mut_N501S;E10_mut_A506V;E10_mut_N501S_A506V];
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

FFU2b = fitness2b./fitness2b(1);
FFU2b_approx = fitness2b_approx./fitness2b_approx(1);


% figure;
% plot([EnergySeq2b-EnergySeq2b(1)],FFU2b,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r2b,p2b] = corr([EnergySeq2b-EnergySeq2b(1)].',FFU2b.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r2b,p2b),'interpreter','tex')

FFU2b_norm = standardize_data(fitness2b);
Energy2b_norm = standardize_data(EnergySeq2b);

FFU2b_approx_norm = standardize_data(fitness2b_approx);
Energy2b_approx_norm = standardize_data(EnergySeq2b_approx);

% figure;
% plot(Energy2b_norm,FFU2b_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r2b,p2b] = corr(Energy2b_norm.',FFU2b_norm.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r2b,p2b),'interpreter','tex')
% 
% figure;
% plot(Energy2b_approx_norm,FFU2b_approx_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r2b_approx,p2b_approx] = corr(Energy2b_approx_norm.',FFU2b_approx_norm.','type','spearman');
% title(sprintf('Spearman r = %.3f, p = %.3e',r2b_approx,p2b_approx),'interpreter','tex')



%% Paper 3: [Pierce2016]

[~,~,dataInput] = xlsread('Infectivity_E2data_Pierce2016.xlsx');

for kk = 1:size(dataInput,1)-1
    strainName3{kk} = dataInput{kk+1,1};
%     if strcmp(strainName{kk},'H77')~=1
%        amino_orig(kk) =  strainName{kk}(1);
%        amino_mut(kk) =  strainName{kk}(end);
%     end
    mut_site3(kk) = str2double(dataInput{kk+1,1}(end-2:end));
    amino_mut3(kk) =  strainName3{kk}(1);
    fitness3(kk) = dataInput{kk+1,2};
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



FFU3_norm = standardize_data(fitness3([1 3 5 11 16 18:20 24 25 34 39 41 43 46 48 62:63]));
%significant (*) and greater than detection threshold
FFU3_norm(8) = [];
EnergySeq3 = EnergySeq3([1 3 5 11 16 18:20 24 25 34 39 41 43 46 48 62:63]); 
EnergySeq3(8) = [];  %mutation A403 doesnt exist 
Energy3_norm = standardize_data(EnergySeq3);

% figure;
% plot(Energy3_norm,FFU3_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r3,p3] = corr(Energy3_norm.',FFU3_norm.','type','spearman');
% title(sprintf('Spearman r = %.3f, p = %.3e',r3,p3),'interpreter','tex')


%% Paper 6 [Goffard2005]

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
H77_mut_N476Q(rev_translation_indices(476-383,ind_non_conserve)) = 'G';

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
H77_mut_N576Q(rev_translation_indices(576-383,ind_non_conserve)) = 'S';

H77_mut_N623Q = H77_mut;
H77_mut_N623Q(rev_translation_indices(623-383,ind_non_conserve));
mutant_order{rev_translation_indices(623-383,ind_non_conserve)};
H77_mut_N623Q(rev_translation_indices(623-383,ind_non_conserve)) = 'D';

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

FFU6 = FFU6([1:7 9]); 
%removing these points as less infectivity was due to lack of incorporation of E1E2 heterodimer into HCVpp.
EnergySeq6 = EnergySeq6([1:7 9]);

% figure;
% plot([EnergySeq6-EnergySeq6(1)],FFU6,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r6,p6] = corr([EnergySeq6-EnergySeq6(1)].',FFU6.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r6,p6),'interpreter','tex')

FFU6_norm = standardize_data(fitness6([1:7 9]));
Energy6_norm = standardize_data(EnergySeq6);

% figure;
% plot(Energy6_norm,FFU6_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r6,p6] = corr(Energy6_norm.',FFU6_norm.','type','spearman');
% title(sprintf('Spearman r = %.3f, p = %.3e',r6,p6),'interpreter','tex')

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

fitness7 = [100 21.68674699	9.036144578	11.44578313	5.421686747	61.44578313	6.024096386	94.57831325	2.409638554	3.012048193	33.73493976	85.54216867	1.204819277	88.55421687	10.84337349 4.81927710843373];

FFU7 = fitness7./fitness7(1);

FFU7_norm = standardize_data(fitness7);
Energy7_norm = standardize_data(EnergySeq7);

% figure;
% plot(Energy7_norm,FFU7_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r7,p7] = corr(Energy7_norm.',FFU7_norm.','type','spearman');
% title(sprintf('Spearman r = %.3f, p = %.3e',r7,p7),'interpreter','tex')

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
EnergySeq8 = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq8(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end

EnergySeq8(9) = [];%(EnergySeq8(1)-min(diag(H)));
fitness8 = [100 77 64 2 2 1 2 4 5 1 12 2 2 19 3 3 44 2 2 2]; %Table 1 Huh7
fitness8_approx = [fitness8(1:3)  mean(fitness8([4:10 12 13 15 16 18:20])) fitness8([11 14 17])];
EnergySeq8_approx = [EnergySeq8(1:3)  mean(EnergySeq8([4:10 12 13 15 16 18:20])) EnergySeq8([11 14 17])];

FFU8 = fitness8./fitness8(1);
FFU8_approx = fitness8_approx./fitness8_approx(1);

% figure;
% plot([EnergySeq8-EnergySeq8(1)],FFU8,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r8,p8] = corr([EnergySeq8-EnergySeq8(1)].',FFU8.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r8,p8),'interpreter','tex')

FFU8_norm = standardize_data(fitness8);
Energy8_norm = standardize_data(EnergySeq8);

% figure;
% plot(Energy8_norm,FFU8_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r8,p8] = corr(Energy8_norm.',FFU8_norm.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r8,p8),'interpreter','tex')

FFU8_approx_norm = standardize_data(fitness8_approx);
Energy8_approx_norm = standardize_data(EnergySeq8_approx);

% figure;
% plot(Energy8_approx_norm,FFU8_approx_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r8_approx,p8_approx] = corr(Energy8_approx_norm.',FFU8_approx_norm.','type','spearman');
% title(sprintf('Spearman r = %.3f, p = %.3e',r8_approx,p8_approx),'interpreter','tex')

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

H77_mut_A439G = H77_mut;
H77_mut_A439G(rev_translation_indices(439-383,ind_non_conserve));
mutant_order{rev_translation_indices(439-383,ind_non_conserve)};
H77_mut_A439G(rev_translation_indices(439-383,ind_non_conserve)) = 'T';

H77_mut_G440A = H77_mut;
H77_mut_G440A(rev_translation_indices(440-383,ind_non_conserve));
mutant_order{rev_translation_indices(440-383,ind_non_conserve)};
H77_mut_G440A(rev_translation_indices(440-383,ind_non_conserve)) = 'A';

H77_mut_L441A = H77_mut;
H77_mut_L441A(rev_translation_indices(441-383,ind_non_conserve));
mutant_order{rev_translation_indices(441-383,ind_non_conserve)};
H77_mut_L441A(rev_translation_indices(441-383,ind_non_conserve)) = 'P';

H77_mut_F442L = H77_mut;
H77_mut_F442L(rev_translation_indices(442-383,ind_non_conserve));
mutant_order{rev_translation_indices(442-383,ind_non_conserve)};
H77_mut_F442L(rev_translation_indices(442-383,ind_non_conserve)) = 'L';

H77_mut_Y443S = H77_mut;
H77_mut_Y443S(rev_translation_indices(443-383,ind_non_conserve));
mutant_order{rev_translation_indices(443-383,ind_non_conserve)};
H77_mut_Y443S(rev_translation_indices(443-383,ind_non_conserve)) = 'S';

msaToTest = [H77_mut;H77_mut_G436A;H77_mut_W437F;H77_mut_L438V;H77_mut_A439G;...
    H77_mut_G440A;H77_mut_L441A;H77_mut_F442L;H77_mut_Y443S];
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

fitness9 = [100 10.31746032	40.21164021	12.6984127	10	19.57671958	14.28571429	62.96296296	6.878306878];

FFU9 = fitness9./fitness9(1);

% figure;
% plot([EnergySeq9-EnergySeq9(1)],FFU9,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r9,p9] = corr([EnergySeq9-EnergySeq9(1)].',FFU9.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r9,p9),'interpreter','tex')

FFU9_norm = standardize_data(fitness9);
Energy9_norm = standardize_data(EnergySeq9);

% figure;
% plot(Energy9_norm,FFU9_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r9,p9] = corr(Energy9_norm.',FFU9_norm.','type','spearman');
% title(sprintf('Spearman r = %.3f, p = %.3e',r9,p9),'interpreter','tex')

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

% figure;
% plot([EnergySeq10-EnergySeq10(1)],FFU10,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r10,p10] = corr([EnergySeq10-EnergySeq10(1)].',FFU10.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r10,p10),'interpreter','tex')

FFU10_norm = standardize_data(fitness10);
Energy10_norm = standardize_data(EnergySeq10);

% figure;
% plot(Energy10_norm,FFU10_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
% [r10,p10] = corr(Energy10_norm.',FFU10_norm.','type','spearman')
% title(sprintf('Spearman r = %.3f, p = %.3e',r10,p10),'interpreter','tex')

FFU10_approx_norm = standardize_data(fitness10_approx);
Energy10_approx_norm = standardize_data(EnergySeq10_approx);

% figure;
% plot(Energy10_approx_norm,FFU10_approx_norm,'o','MarkerFaceColor',blue,'MarkerEdgeColor',blue)
[r10_approx,p10_approx] = corr(Energy10_approx_norm.',FFU10_approx_norm.','type','spearman');
% title(sprintf('Spearman r = %.3f, p = %.3e',r10_approx,p10_approx),'interpreter','tex')

%% Fitness vs Energy (normalized) including 02.E10 [merging similar ones]

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

markersize = 9;
line_width = 0.5;

pink = color_scheme(8,:);
figure;
plot(Energy6_norm,FFU6_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',red,'LineWidth',line_width)
% xlabel('E_{mutant} - E_{H77}','interpreter','tex')
% ylabel('f_{mutant} / f_{H77}','interpreter','tex')
xlabel('Energy (normalized)')
ylabel('Experimental fitness (normalized)')
% xlabel('$${\rm E}_{\rm mutant}-{\rm E}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$${\rm f}_{\rm mutant}/{\rm f}_{\rm Mahoney}$$','interpreter','latex')
% ylabel('$$\frac{f_{mutant}}{f_{Mahoney}}$$','interpreter','latex')
hold on 
plot(Energy9_norm,FFU9_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',blue,'LineWidth',line_width)
plot(Energy7_norm,FFU7_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',green,'LineWidth',line_width)
plot(Energy8_approx_norm,FFU8_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'LineWidth',line_width)
% plot(Energy8b,FFU8b,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',purple,'LineWidth',line_width)
plot(Energy1_norm,FFU1_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',orange,'Linewidth',line_width)
% plot(Energy5,titer5,'<','MarkerEdgeColor',purple,'MarkerSize',10,'Linewidth',1.5)
plot(Energy2_norm,FFU2_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',yellow,'LineWidth',line_width)
plot(Energy2b_approx_norm,FFU2b_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',brown,'Linewidth',line_width)%,'MarkerFaceColor','r')
plot(Energy10_approx_norm,FFU10_approx_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',pink,'Linewidth',line_width)
plot(Energy3_norm,FFU3_norm,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor','k','Linewidth',line_width)
% plot(Energy2,titer2,'o','MarkerSize',markersize,'MarkerEdgeColor','w','MarkerFaceColor',lightgreen,'Linewidth',line_width)
% legend('Gal-Tanamy2008','Keck2009-H77','Keck2009-02.E10','Pierce2016')
h = legend('Goffard2005','Drummer2006','Falkowska2007','Rothwangl2008','Gal-Tanamy2008','Keck2009','Keck2009b','Keck2012','Pierce2016','Location','NorthEast');
h.Position = [0.70 0.75 0.12 0.0488];
legend boxoff
ylim([-2.5 4.5])
xlim([-3.5 3.5])
% xlabel('E-Eref')
% ylabel('f/fref');

% [rho,pval] = corr([Energy1 Energy2 Energy2b Energy3]',[FFU1 FFU2 FFU2b FFU3]','tail','left','type','spearman')
% [rho,pval] = corr([Energy10_approx_norm Energy9_norm Energy8_approx_norm Energy7_norm Energy6_norm Energy1_norm Energy2_norm Energy2b_approx_norm Energy3_norm]',...
%     [FFU10_approx_norm FFU9_norm FFU8_approx_norm FFU7_norm FFU6_norm FFU1_norm FFU2_norm FFU2b_approx_norm FFU3_norm]','tail','left','type','spearman');

% rho_averge = mean([r10_approx r9 r8_approx r7 r6 r1 r2 r3 r2b_approx])
% pval_average = mean([p10_approx p9 p8_approx p7 p6 p1 p2 p3 p2b_approx])

total_length = (length(FFU2b_approx_norm)+length(FFU6_norm)+length(FFU1_norm)+length(FFU2_norm)+...
    length(FFU3_norm)+length(FFU7_norm)+length(FFU8_approx_norm)+length(FFU9_norm)+length(FFU10_approx_norm));
w1 = length(FFU1_norm)/total_length;
w2 = length(FFU2_norm)/total_length;
w3 = length(FFU3_norm)/total_length;
w6 = length(FFU6_norm)/total_length;
w7 = length(FFU7_norm)/total_length;
w8 = length(FFU8_approx_norm)/total_length;
w9 = length(FFU9_norm)/total_length;
w2b = length(FFU2b_approx_norm)/total_length;
w10 = length(FFU10_approx_norm)/total_length;


rho_weighted_average = r1*w1+r2*w2+r3*w3+r6*w6+r7*w7+r8_approx*w8+r9*w9+r2b_approx*w2b+r10_approx*w10;

% title(sprintf('r = %0.3f, p = %0.3e',rho,pval))
text(-3,-2,sprintf('$$\\bar{r} = %.2f$$',rho_weighted_average),'interpreter','latex','FontSize',16)

% P = polyfit([Energy1 Energy2 Energy2b Energy3]',[FFU1 FFU2 FFU2b FFU3]',1);
P = polyfit([Energy1_norm Energy2_norm Energy2b_approx_norm Energy3_norm Energy6_norm Energy7_norm Energy8_approx_norm Energy9_norm Energy10_approx_norm]',...
    [FFU1_norm FFU2_norm FFU2b_approx_norm FFU3_norm FFU6_norm FFU7_norm FFU8_approx_norm FFU9_norm FFU10_approx_norm]',1);
x = -3:.5:3; %xaxis
y = P(1)*x+P(2);
plot(x,y,'k--','LineWidth',1)

set(gca, ...
  'Box'         , 'off'     , ...
  'TickDir'     , 'out'     , ...
  'TickLength'  , [.02 .02] , ...
  'XMinorTick'  , 'off'      , ...
  'YMinorTick'  , 'off'      , ...
  'XColor'      , [.2 .2 .2], ...
  'YColor'      , [.2 .2 .2], ...
  'YTick'       , -3:1:4, ...
  'XTick'       , -3:1:3, ...
  'LineWidth'   , 0.5        );
