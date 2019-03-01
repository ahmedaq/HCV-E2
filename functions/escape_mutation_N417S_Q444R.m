function escape_mutation_N417S_Q444R(msa_aa, phi_curr, mutant_order, ...
    H, ind_non_conserve)

% Code for analyzing HCV1 escape mutation N417S and associated compensatory
% mutation Q444R
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-08-24

%%

run startup.m

H = convertRayFLParamsToJohnParamsFormat(H);

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

%% Morin 2012 extended study of compensatory mutation available for N417S mutation

set(0,'DefaultTextFontSize',8)
set(0,'DefaultAxesFontSize',8)

H77_mut_N417S = H77_mut;
H77_mut_N417S(rev_translation_indices(417-383,ind_non_conserve));
mutant_order{rev_translation_indices(417-383,ind_non_conserve)};
H77_mut_N417S(rev_translation_indices(417-383,ind_non_conserve)) = 'S';
baseSeq = H77_mut_N417S;

muts_Babcock = [603 563 624 394 481 408 618 399 434 463 558];
indx_muts_Babcock = [];

nn=0;
for kk = 1:length(mutant_order)
    if kk~=417-383
        muts_possible{kk} = setdiff(mutant_order{kk},baseSeq(kk));
        for mm = 1:length(muts_possible{kk})
            nn = nn+1;
            seq_test = baseSeq;
            seq_test(kk) = muts_possible{kk}(mm);
            seqs_with_mutants_in_baseSeq(nn,:) = seq_test;
            mut_identity{nn} = sprintf('%s%d%s',baseSeq(kk),ind_non_conserve(kk)+383,seq_test(kk));
            if ismember(ind_non_conserve(kk)+383, muts_Babcock)
                indx_muts_Babcock = [indx_muts_Babcock nn];
            end
        end
    end
end

msaToTest = [baseSeq;H77_mut;seqs_with_mutants_in_baseSeq];
numSeqMsaToTest = size(msaToTest,1);
EnergySeq_N417S = zeros(1,numSeqMsaToTest);
for i = 1:numSeqMsaToTest
    input_seq = msaToTest(i,:);
    
    input_parm = cell(1,2);
    input_parm{1} = size(msaToTest,2);
    input_parm{2} = input_seq;
    
    [out_seq_ex(i,:)] = convertAAseq2Bin_new(mutant_order,bin_matrix,input_parm);
    [EnergySeq_N417S(i)] = calcSeqEnergy(out_seq_ex(i,:),H);
end     

Energy_baseSeq = EnergySeq_N417S(1);
% Energy_H77mut = EnergySeq_N417S(2);
EnergySeq_N417S(1:2) = [];
        
% mut_identity(find(EnergySeq_N417S-Energy_baseSeq<0)) 
% find(EnergySeq_N417S-Energy_baseSeq<0)
% indx_Q444R = find(strcmp(mut_identity,'Q444R'))
% EnergySeq_N417S(find(EnergySeq_N417S-Energy_baseSeq<0))-Energy_baseSeq

figure;
gray = [0.6 0.6 0.6];
generatePDF([EnergySeq_N417S-Energy_baseSeq].',gray,'count',50)
set(gca,'XMinorTick','off',...
    'YMinorTick','off')
xlabel('Energy (H77_{N417S+X}) - Energy (H77_{N417S})','FontSize',8)
xlim([-6 12])
