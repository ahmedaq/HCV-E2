function JohnFLParamsMatrix = convertRayFLParamsToJohnParamsFormat(RayFLParamsMatrix)

% Code for converting the model parameters matrix inferred using MPF-BML
% method into the format of model parameters matrix generated using ACE
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%

diag_RayFLParamsMatrix = diag(diag(RayFLParamsMatrix));
offdiag_RayFLParamsMatrix = RayFLParamsMatrix - diag_RayFLParamsMatrix;
offdiag_JohnFLParamsMatrix = 2*offdiag_RayFLParamsMatrix;
JohnFLParamsMatrix = -(offdiag_JohnFLParamsMatrix + diag_RayFLParamsMatrix);

for kk = 1:length(JohnFLParamsMatrix)
    JohnFLParamsMatrix(kk+1:end,kk) = 0;
end