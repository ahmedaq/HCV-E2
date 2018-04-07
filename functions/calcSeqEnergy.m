function E = calcSeqEnergy(seq,H)

% Code for computing the sequence energy using the model parameters
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

E = seq*(-triu(H))*seq';