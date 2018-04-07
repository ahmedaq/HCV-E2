function plot_triple_mutant_probability(p3_data,p3_sampler)

% Code for plotting the triple mutant probabilities
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%

run startup.m
set(0,'DefaultTextFontSize',10)
set(0,'DefaultAxesFontSize',10)

color = [0.6000    0.6000    0.6000];%gray
markersize = 3;

figure
plot(p3_data,p3_sampler,'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize);hold on;
max_value = max(max(p3_data,p3_sampler));
plot(0:max_value/20:max_value,0:max_value/20:max_value,'k')
xlabel('Triple mutant probability (MSA)')
ylabel('Triple mutant probability')
[r,pval] = corr(p3_data.',p3_sampler.','type','pearson','tail','both');
% title( ['\epsilon_c = ' num2str(sum_error_nondiag) ] )
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'XTick'       , 0:0.1:0.3, ...
    'YTick'       , 0:0.1:0.3, ...
    'LineWidth'   , .5       );