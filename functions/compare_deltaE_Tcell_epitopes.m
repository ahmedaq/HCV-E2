function compare_deltaE_Tcell_epitopes(dE2,ind_non_conserve)

% Code for comparing fitness costs associated with residues defining known
% T cell epitopes
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%
run startup.m

%% NW9
gray = [0.6 0.6 0.6];
NW9 = [541:550];
figure_epitopes2(-dE2,NW9,ind_non_conserve,'',gray,'Epitope associated with spontaneous clearance',0)
ylim([0 14])
post_proc_fig_box

%% Analysis of all CD8+ epitopes 

clear CD8_ep

CD8_ep{1} = [489:496]; %B51
CD8_ep{2} = [530:539]; %B60
CD8_ep{3} = [541:550]; %B57 clearance
CD8_ep{4} = [632:641]; %A3 clearance
CD8_ep{5} = [630:639]; %A3 clearance

CD8_ep{6} = [401:411]; %A2
CD8_ep{7} = [459:469]; %B53

CD8_ep{8} = [497:507]; %B35 [Ward2002]

CD8_ep{9} = [569:578]; %B50
CD8_ep{10} = [610:618]; %cw7
CD8_ep{11} = [613:622]; %cw7, A2
CD8_ep{12} = [621:628]; %A11 clearance
CD8_ep{13} = [654:662]; %B60
CD8_ep{14} = [686:694]; %A2
CD8_ep{15} = [712:726]; %A*2402
CD8_ep{16} = [723:734]; %A2


figure;
G = [];
data = [];
data_cons = [];
data_entropy = [];

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)


% boxplot

eps_control = [CD8_ep{3:5}];
eps_rest = [CD8_ep{[1:2 6:16]}];

G = [zeros(1,length(eps_control)) ...
    ones(1,length(eps_rest))];
data = [dE2(eps_control-383) ...
    dE2(eps_rest-383)];

box_lineWidth = 0.5;
box_widths_value = 0.2;
box_color = [0.6 0.6 0.6];
box_color_transparency = 0.4; %faceAlpha
median_lineWidth = 1;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = 'o';
outlier_markerSize = 4;
outlier_marker_edgeWidth = 0.1;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0.75;
label_xaxis_data = {'Spontaneous clearance','Remaining'};
text_ylabel = '\DeltaE_i';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = -1;
ylim_max = 14;
savefig = 0;
savefig_name = 'CTL_control_rest';
fig_width_cm = 4;
fig_height_cm = 5;

figure
figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);


p_diff_dist_control_rest = ranksum(dE2(eps_control-383),dE2(eps_rest-383));
fprintf('\nP = %.1e, Mann-Whitney test\n',p_diff_dist_control_rest)
