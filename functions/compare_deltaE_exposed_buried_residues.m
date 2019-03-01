function compare_deltaE_exposed_buried_residues(dE2)

% Code for comparing fitness costs associated with mutations in the 
% exposed and buried residues
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%
area_each_atom = xlsread('area_each_atom_4mwf.xlsx'); %using get_area() in Pymol (PDB ID: 4MWF)
load order_E2_4mwf.mat
order_unique = unique(order);

%Accessible surface area (ASA) of residue = sum of accesible surface area of each
%atom
for kk = 1:length(order_unique)
    ASA(kk) = sum(area_each_atom(find(order == order_unique(kk)))); %residue level
end

%Finding MaxASA for each residue in 4MWF
load res_nos_E2_4mwf.mat 
load H77_seq
seq_4mwf = H77(unique(order)-383);
%two changes in structure in pymol from H77
seq_4mwf(find(order_unique==448)) = 'D'; 
seq_4mwf(find(order_unique==576)) = 'D';

[A,B,C] = xlsread('MaxASA_residues.xlsx'); %from https://en.wikipedia.org/wiki/Relative_accessible_surface_area
MaxASA_Miller1987 = A(:,3);
aminoacids = [B{:,2}];

for kk = 1:length(seq_4mwf)
    MaxASA_seq_4mwf_Miller1987(kk)=MaxASA_Miller1987(find(ismember(aminoacids,seq_4mwf(kk))));
end
    
RASA_4mwf_Miller1987 = ASA./MaxASA_seq_4mwf_Miller1987;

thresh_RASA = 0.2; 

buried_residues = res_nos2(find(RASA_4mwf_Miller1987<=thresh_RASA));

rest = [384:410 setdiff(res_nos2,buried_residues)]; %including HVR1 too

dE2_4mwf_buried_residues_get_area = dE2(buried_residues - 383);
dE2_4mwf_surface_residues_get_area = dE2(rest - 383);

p_mannwhitney_surface_buried = ranksum(dE2_4mwf_buried_residues_get_area,dE2_4mwf_surface_residues_get_area,'tail','right');
fprintf('\nP = %.1e, Mann-Whitney test\n',p_mannwhitney_surface_buried)


% BOXPLOT
G = [];
data = [];

figure
G = [zeros(1,length(dE2_4mwf_buried_residues_get_area)) ...
    ones(1,length(dE2_4mwf_surface_residues_get_area))];
data = [(dE2_4mwf_buried_residues_get_area) ...
    (dE2_4mwf_surface_residues_get_area)];

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

box_lineWidth = 0.5;
box_widths_value = 0.2;
box_color = gray;
box_color_transparency = 0.2; %faceAlpha
median_lineWidth = 1;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = 'o';
outlier_markerSize = 8;
outlier_marker_edgeWidth = 0.1;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;
label_xaxis_data = {'Buried','Exposed'};
text_ylabel = 'Escape time';
text_xlabel = '';
text_title = '';%'E2 surface and buried residues';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 1000;
savefig = 0;
savefig_name = 'surface_buried_residues';
fig_width_cm = 4;
fig_height_cm = 5;

figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);