function compare_deltaE_known_escape_mutations(dE2,L)

% Code for comparing fitness costs associated with known escape mutations
% and those on remaining residues in E2
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07


%%
run startup.m

%Escape mutations [Bailey 2015]
polymorphisms_associated_with_neutralization_resistance1 = ... 
    [416 422 424 431 433 438 442 446 453 456 461 475 482 520 524 531 533 557 558 560];

%Escape mutations with reasonable infectivity [Keck2009]
polymorphisms_associated_with_neutralization_resistance2 = ...
    [431 435 466 528 531 538 580 610 636 713 444 446 482 501]; %type 1 (fixed in 1991), see Keck2009

%Escape mutations [Morin2012]
polymorphisms_associated_with_neutralization_resistance3 = ...
    [415 417 444];

%Escape mutation [Keck2008]
polymorphisms_associated_with_neutralization_resistance4 = 431;

polymorphisms_associated_with_neutralization_resistance = unique([...
    polymorphisms_associated_with_neutralization_resistance1 ...
    polymorphisms_associated_with_neutralization_resistance2 ...
    polymorphisms_associated_with_neutralization_resistance3 ...
    polymorphisms_associated_with_neutralization_resistance4]);

%boxplot

G = [zeros(1,length(polymorphisms_associated_with_neutralization_resistance)) ...
    ones(1,length(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)))];
data = [dE2(polymorphisms_associated_with_neutralization_resistance-383) ...
    dE2(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383))];

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

box_lineWidth = 0.5;
box_widths_value = 0.2;
box_color = [0 0 0];
box_color_transparency = 0.2; %faceAlpha
median_lineWidth = 1;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = 'o';
outlier_markerSize = 4;
outlier_marker_edgeWidth = 0.1;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0.75;
label_xaxis_data = {'Escape','Remaining'};
text_ylabel = 'Fitness Costs, \DeltaE_i';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = -1;
ylim_max = 14;
savefig = 0;
savefig_name = 'escape_mutations';
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

P = ranksum(dE2(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),dE2(polymorphisms_associated_with_neutralization_resistance-383));
fprintf('\nP = %.1e, Mann-Whitney test\n',P)
% [h,p] = kstest2(dE2(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),dE2(polymorphisms_associated_with_neutralization_resistance-383),'tail','smaller')
