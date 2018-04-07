
function [freqCountNumMutPerSeq numSeq] = plotXMutFreqPMF(mutMtxIn, protLen)

% Code for finding frequency of each mutant
% 
% Written by: Muhammad Saqib Sohail
% Edited by: Ahmed Abdul Quadeer
% Last updated: 2018-04-07

%%

numMutPerSeq = sum(mutMtxIn, 2);

numSeq = size(mutMtxIn,1);

freqCountNumMutPerSeq = zeros(1,protLen);

for i = 1:protLen+1
   freqCountNumMutPerSeq(i) = sum(numMutPerSeq == i-1) ;
    
end

return;