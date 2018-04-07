
function [xMutMtxCell, allMutMtx] = mutCountMSA(msa_aa_ex, phi_curr, protLen)

% Code for finding number of mutations per sequence
% 
% Written by: Muhammad Saqib Sohail
% Edited by: Ahmed Abdul Quadeer
% Last updated: 2018-04-07

%%
numSeq = size(msa_aa_ex, 1);
msa_aa_exLen = size(msa_aa_ex, 2);

maxNumMutPerSite = max(phi_curr);

% each cell location witll contain a numSeq x protLen matx of 1 and 0 with
% 1 indicating a mutation. the xth cell entry will contain the mut msa for
% the xth dominant mutant.
xMutMtxCell = cell(1, maxNumMutPerSite); 

zeroCol = zeros(numSeq, 1);

for i = 1:maxNumMutPerSite
%     i
%     phi_curr(1:5)
    lastCol = 0;
    mutMtx_temp = zeros(numSeq, protLen);
    for j = 1:protLen
%         j
        numCol = phi_curr(j);
        actualColNumber = numCol + 1 - i;
                
        if(actualColNumber > 0)
            mutMtx_temp(:,j) = msa_aa_ex(:, actualColNumber + lastCol);
            
%                 if(j == 275)
%                     display('---------------')
%                     i
%                     j
%                     numCol
%                     actualColNumber
%                     lastCol
%                     co = actualColNumber + lastCol
%                     msa_aa_ex(896, actualColNumber + lastCol)
%                     display('RUNNING IF')
%                     display('---------------')
%                 end
            
        else
            
%                 if(j == 275)
%                     display('++++++++++++++++++++++')
%                     i
%                     j
%                     numCol
%                     actualColNumber
%                     %msa_aa_ex(896, actualColNumber + lastCol)
%                     display('RUNNING ELSE')
%                     display('++++++++++++++++++++++')
%                 end
            
            mutMtx_temp(:,j) = zeroCol;
        end
        lastCol = sum(phi_curr(1:j));
%         pause
    end
    xMutMtxCell{i} = mutMtx_temp;
end

allMutMtx = zeros(numSeq, protLen);
for i = 1:maxNumMutPerSite
    allMutMtx = allMutMtx + xMutMtxCell{i};
    if(max(max(allMutMtx)) >1 )
        i
        pause
    end
end
