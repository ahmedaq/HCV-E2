function mapping_crystal_structure_pymol(mean_escape_time)

% Code for mapping fitness costs on E2 crystal structure
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-12-21

load res_nos_E2_4mwf.mat %file containing the positions of E2 for which structure is available. obtained from 4mwf structure code in polio folder
load order_E2_4mwf.mat %file containing the numbering of all atoms (not just carbon alpha)


%% making a vector of length(order) for assigning dE values to all atoms of residues

for kk = 1:length(order)
    escape_time_crystal_structure_4mwf_allatoms(kk) = mean_escape_time(order(kk)-383);
end
delete crystal_structure_4mwf_allatoms_mean_escape_time.txt
dlmwrite('crystal_structure_4mwf_allatoms_mean_escape_time.txt',escape_time_crystal_structure_4mwf_allatoms.')

%% Defining regions for pymol 

hvr1 = 384:408;
cd81_binding = [420 421 424 427 430 432 436:438 440:443 523 526 527 529 530 535 540 549 550 613 614 616:618]; %Pierce2016 [TableS5]
%CBH-23
HmAb_Pierce2016{10} = [494 508 509 537 539 549 552 554 564 611 614 644]; %same as CBH7
%HC33-1
HmAb_Pierce2016{15} = [413 418 420];

pymol_input_hvr1 = pymol_data_input(hvr1,'HVR1');
pymol_input_cd81 = pymol_data_input(cd81_binding,'cd81bs');
pymol_input_HC331 = pymol_data_input(HmAb_Pierce2016{15},'HC331');
pymol_input_CBH23 = pymol_data_input(HmAb_Pierce2016{10},'CBH23');

pymol_input_hvr1
pymol_input_cd81
pymol_input_HC331
pymol_input_CBH23


