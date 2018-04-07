function msa_new = construct_msa_aa_after_entropy_compression(msa,mutant_order)

% Code for constructing a modified MSA based on the entropy reduction done
% in MPF-BML
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

% Inputs:
% msa = original amino acid msa without entropy compression
% mutant_order = list of amino acids at each site after entropy compression
% Output:
% msa_new = modified amino acid msa with only those amino acids at each
%           site that are kept in mutant_order (after entropy compression)

%%

[ns,ls] = size(msa);

msa_new = msa;
for kk = 1:ls
    indices_aa_in_mutant_order = [];
    for mm = 1:length(mutant_order{kk})
        indices_aa_in_mutant_order = [indices_aa_in_mutant_order;find(msa(:,kk)==mutant_order{kk}(mm))];
    end
    indices_aa_not_in_mutant_order = setdiff(1:ns,indices_aa_in_mutant_order);
    msa_new(indices_aa_not_in_mutant_order,kk) = mutant_order{kk}(end);
end
