function Out_sector = rev_translation_indices(In_sector,V)

% Code for translating the indexing from actual (length L protein) to only
% mutant number of resdiues protein (length L_mut)
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%

Out_sector = [];

for mm = 1:length(In_sector)
    for kk = 1:length(V)
        if V(kk)==In_sector(mm)
            Out_sector(mm)=kk;
        end
    end
end
    


