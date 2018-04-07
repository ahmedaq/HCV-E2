function pymol_input = pymol_data_input(sites_vector,name)

% Code for generating code to be used in Pymol for selecting a group of
% residues
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

pymol_input_temp = sprintf('+%d',sites_vector);
pymol_input = sprintf('select %s, resi %s',name,pymol_input_temp);