function compare_deltaE_antigenic_domains(mean_escape_time)

% Code for comparing fitness costs associated with regions targeted by
% HmAbs
%
% Written by: Ahmed Abdul Quadeer
% Last updated: 2018-12-21

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

%% Bar plot -- minimume escape time

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

for kk = 1:length(antigenic_region)
    data_mean{kk} = mean_escape_time(antigenic_region{kk}-383);
    min_data_mean(kk) = min(data_mean{kk});
end


ylim_min = 0;
ylim_max = 300;

figure
subplot(2,1,1)
hold on
xbars = [6.5 7.5];
patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+1 ylim_max ylim_max ylim_min+1], [0.85 0.85 0.85], ...
    'EdgeColor','w','LineWidth',0.1)

xbars = [4.5 5.5];
patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+1 ylim_max ylim_max ylim_min+1], [0.85 0.85 0.85], ...
    'EdgeColor','w','LineWidth',0.1)

bar(1:7, [zeros(1,0) min_data_mean(1) zeros(1,6)], 0.4, 'FaceColor','k','FaceAlpha',0.7);
bar(1:7, [zeros(1,1) min_data_mean(2) zeros(1,5)], 0.4, 'FaceColor',darkgray,'FaceAlpha',0.7);
bar(1:7, [zeros(1,2) min_data_mean(3) zeros(1,4)], 0.4, 'FaceColor',color_scheme_npg(1,:),'FaceAlpha',0.7);
bar(1:7, [zeros(1,3) min_data_mean(4) zeros(1,3)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);
bar(1:7, [zeros(1,4) min_data_mean(5) zeros(1,2)], 0.4, 'FaceColor',color_scheme_npg(3,:),'FaceAlpha',0.7);
bar(1:7, [zeros(1,5) min_data_mean(6) zeros(1,1)], 0.4, 'FaceColor',color_scheme_npg(4,:),'FaceAlpha',0.7);
bar(1:7, [zeros(1,6) min_data_mean(7) zeros(1,0)], 0.4, 'FaceColor',color_scheme_npg(5,:),'FaceAlpha',0.7);


set(gca,'XTickLabel',{'HVR1','CD81bs','Domain A','Domain B','Domain C','Domain D','Domain E'})
xtickangle(30)
ylabel('Minimum escape time')
xlim([0.5 7.5])
post_proc_fig

%% Heatmap -- for calculating n_e^tau

thresh_range = 80:20:140;
for mm = 1:length(thresh_range)
    for kk = 1:length(data_mean)
        no_of_sites_less_100(kk,mm) = sum(data_mean{kk}<=thresh_range(mm));
    end
end

labels_cell{2} = {'100','140','180','220'};
labels_cell{1} = {'HVR1','CD81bs','Domain A','Domain B','Domain C','Domain D','Domain E'};

% figure_heatmap(no_of_sites_less_100_HmAb20.','BuPu','',{'','Escape time'},[0 8],labels_cell)

subplot(2,1,2)

heatmap(labels_cell{1},labels_cell{2},no_of_sites_less_100.');
set(gca,'Colormap',color_scheme_BuPu,'ColorLimits',[0 22],'GridVisible','off','FontSize',10,'ColorbarVisible','off');
% h = set(gca,'Colormap',color_scheme,'ColorLimits',limits_data);
% xlabel(text_labels{1})
ylabel('Escape time threshold, tau')
% colorbar('Location','north')
