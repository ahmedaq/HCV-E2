function compare_deltaE_antigenic_domains(dE2)

% Code for comparing fitness costs associated with regions targeted by
% HmAbs
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%
run startup.m

cd81_binding_Castelli2017 = [412:423 424:453 519:535];
cd81_binding_Kong2016 = [420 436:443 527 529 530 535];
cd81_binding_Pierce2016 = [420 421 424 427 430 432 436:438 440:443 523 526 527 529 530 535 540 549 550 613 614 616:618];

antigenic_region = [];

antigenic_region{1} = 384:410; %hvr1
antigenic_region{2} = cd81_binding_Pierce2016; %cd81bs

%based on [Pierce2016], Table 2 caption [previous reported epitopes]
antigenic_region{3} = [626:632]; %domainA
antigenic_region{4} = [425:428 436:438 529 530 535];%%domainB
antigenic_region{5} = [542 544 545 549 561 592 593 598 631 633]; %domainC 
antigenic_region{6} = [441:443 502 616]; %domainD
antigenic_region{7} = 412:423; %domainE

%
G = [];
data = [];

for kk = 1:length(antigenic_region)
    G = [G kk*ones(1,length(antigenic_region{kk}))];    
    data = [data dE2(antigenic_region{kk}-383)];
end

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

black = [0 0 0];

box_lineWidth = 0.5;
box_widths_value = 0.2;
box_color = [black; darkgray; color_scheme_npg(1,:); color_scheme_npg(2,:); ...
    color_scheme_npg(3,:); color_scheme_npg(4,:); color_scheme_npg(5,:)];
% box_color = [black; gray; green; purple; orange; yellow; brown];
box_color_transparency = 0.6; %faceAlpha
median_lineWidth = 1;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = 'o';
outlier_markerSize = 7;
outlier_marker_edgeWidth = 0.1;
outlier_marker_edgeColor = 'k';
outlier_jitter_value = 0.5;
label_xaxis_data = {'HVR1','CD81bs','Domain A','Domain B','Domain C','Domain D','Domain E'};
text_ylabel = 'Fitness cost, \DeltaE_i';
text_xlabel = '';%Regions of interest in E2';
text_title = '';%'Important E2 regions';
label_orientation_choice = 'horizontal'; %'inline'
ylim_min = -1;
ylim_max = 13;
savefig = 0;
savefig_name = 'important_regions_v2';
fig_width_cm = 9;
fig_height_cm = 3;

figure
xbars = [6.5 7.5];
patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
    'EdgeColor','w','LineWidth',0.1)
hold on
xbars = [4.5 5.5];
patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
    'EdgeColor','w','LineWidth',0.1)

figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);