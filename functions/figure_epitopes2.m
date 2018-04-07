function figure_epitopes2(dE,sites,ind_non_conserve,heading,color,xlabel_data,savefigs,newfig)

% Code for generating bar plots
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%
if nargin ~= 8
    newfig = 1;
end

sites_mut = rev_translation_indices(sites-383,ind_non_conserve);
conserved_pos_epitope = find(sites_mut==0);
sites_mut(sites_mut==0) = [];


set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

if newfig == 1
    figure;
end
b = bar(1:length(sites_mut),-dE(sites_mut),0.5,'FaceColor',color,'EdgeColor','k','FaceAlpha',0.5);
xlim([0.5 length(sites_mut)+.5])
ylim([0 12])
xlabel(xlabel_data)
ylabel('\DeltaE_i')
title(heading)
set(gca, ...
    'XTick', 1:length(sites_mut),...
    'XTickLabel', ind_non_conserve(sites_mut)+383)
xtickangle(45)

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.01 .01] , ...
    'XColor'      , [.2 .2 .2], ...
    'YColor'      , [.2 .2 .2], ...
    'LineWidth'   , 1        );

if savefigs
    if length(sites_mut)<20
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 18 6])
        print(sprintf('fig_bar_epitope_%s.png',heading),'-dpng','-r300')
    else
        set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 40 6])
        print(sprintf('fig_bar_epitope_%s.png',heading),'-dpng','-r300')
    end
end

