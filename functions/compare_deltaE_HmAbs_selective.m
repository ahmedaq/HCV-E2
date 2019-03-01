function compare_deltaE_HmAbs_selective(mean_escape_time)

% Code for comparing fitness costs associated with mutating binding
% residues of HmAbs
%
% Written by: Ahmed Abdul Quadeer
% Last updated: 2018-04-07

run startup.m

set(0,'DefaultAxesFontName','Arial')
set(0,'DefaultTextFontName','Arial')
set(0,'DefaultAxesFontSize',10)
set(0,'DefaultTextFontSize',10)

%% Pierce2016, Table S5, without non-neutralizing AR1A and AR1B v2 (without those in Table S1)
%% Also including other HmAbs
% [Owsianka2008]
thresh_binding_imp = 25; %<%
[dataInput1] = xlsread('Binding_CHB5_Owsianka2008.xlsx');
CBH5 =  dataInput1(dataInput1(:,2)<thresh_binding_imp,1);
CBH5 =  unique([CBH5.' 494 497 614 617:619 621 624 502:505 507:509]); %extra sites in [lacob2008]
[dataInput2] = xlsread('Binding_CHB7_Owsianka2008.xlsx');
CBH7 =  dataInput2(dataInput2(:,2)<thresh_binding_imp,1);
CBH7 = unique([CBH7.' 494 497 614 617:619 621 624]); %extra sites in [lacob2008]

HmAb_Pierce2016_tableS5 = [];

HmAb_Pierce2016_tableS5{1} = [424 523 525 530 535 538 540];%AR3A
HmAb_Pierce2016_tableS5{2} = [412, 416, 418, 423, 424, 523, 525, 530, 535, 540]; %[424 530 535]; %AR3B
HmAb_Pierce2016_tableS5{3} = [424, 488, 523, 525, 530, 535, 538, 540];%[424 525 530 535 540]; %AR3C
HmAb_Pierce2016_tableS5{4} = [412, 424, 523, 530, 535];%[424 530]; %AR3D
HmAb_Pierce2016_tableS5{5} = [412:423]; %HCV1
% HmAb_Pierce2016_tableS5{6} = [413 418 420]; %HC33.1
HmAb_Pierce2016_tableS5{6} = [413 418 420]; %HC33.32
HmAb_Pierce2016_tableS5{7} = [525 530 535]; %HC-1

% HmAb_Pierce2016_tableS5{10} = CBH7; %[Owsianka2008,lacob2008]
HmAb_Pierce2016_tableS5{8} = [523 526 527 529 530 535]; %1:7
HmAb_Pierce2016_tableS5{9} = [523 526 527 529 530 535]; %A8
HmAb_Pierce2016_tableS5{10} = [416 420 529 530 535]; %e137 [Sabo2011]

HmAb_Pierce2016_tableS5{11} = CBH5; %[Owsianka2008,lacob2008]
HmAb_Pierce2016_tableS5{12} = [441 442]; %HC84.1
% HmAb_Pierce2016_tableS5{15} = [441:443 613 616]; %HC84.20
HmAb_Pierce2016_tableS5{13} = [441:443]; %HC84.21
HmAb_Pierce2016_tableS5{14} = [420 428 437 441:443 616]; %HC84.22
HmAb_Pierce2016_tableS5{15} = [420 428 437 441:443 616]; %HC84.23
% HmAb_Pierce2016_tableS5{19} = [442:443]; %HC84.24
HmAb_Pierce2016_tableS5{16} = [441 442 616]; %HC84.25
% HmAb_Pierce2016_tableS5{21} = [441 442 616]; %HC84.26
HmAb_Pierce2016_tableS5{17} = [441:443 446 616]; %HC84.27
% HmAb_Pierce2016_tableS5{23} = [408 413 420]; %HC33.4
HmAb_Pierce2016_tableS5{18} = [408 413 418 420]; %HC33.8
HmAb_Pierce2016_tableS5{19} = [408 413 418 420]; %HC33.29
% HmAb_Pierce2016_tableS5{26} = [425 428 436:438 442 443 530 535]; %HC-11
HmAb_Pierce2016_tableS5{20} = [437 439 530 535 523 431]; %added residues from [Keck2008] %[437 439 530 535]; %CBH-2

data_mean_HmAb_tableS5 = cell(1,length(HmAb_Pierce2016_tableS5));

for kk = 1:length(HmAb_Pierce2016_tableS5)
    data_mean_HmAb_tableS5{kk} = mean_escape_time(HmAb_Pierce2016_tableS5{kk}-383);
    min_data_mean_HmAb_tableS5(kk) = min(data_mean_HmAb_tableS5{kk});
end

ylim_min = 0;
ylim_max = 300;

figure
subplot(2,1,1)

ylim_min = 0;
ylim_max = 300;

xbars = [0.525 10.5];
patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
    'EdgeColor','w','LineWidth',0.1)
hold on

bar(1:20, [zeros(1,0) min_data_mean_HmAb_tableS5(1) zeros(1,19)], 0.4, 'FaceColor',darkgray,'FaceAlpha',0.7);
bar(1:20, [zeros(1,1) min_data_mean_HmAb_tableS5(2) zeros(1,18)], 0.4, 'FaceColor',darkgray,'FaceAlpha',0.7);
bar(1:20, [zeros(1,2) min_data_mean_HmAb_tableS5(3) zeros(1,17)], 0.4, 'FaceColor',darkgray,'FaceAlpha',0.7);
bar(1:20, [zeros(1,3) min_data_mean_HmAb_tableS5(4) zeros(1,16)], 0.4, 'FaceColor',darkgray,'FaceAlpha',0.7);
bar(1:20, [zeros(1,4) min_data_mean_HmAb_tableS5(5) zeros(1,15)], 0.4, 'FaceColor',darkgray,'FaceAlpha',0.7);
bar(1:20, [zeros(1,5) min_data_mean_HmAb_tableS5(6) zeros(1,14)], 0.4, 'FaceColor',darkgray,'FaceAlpha',0.7);
bar(1:20, [zeros(1,6) min_data_mean_HmAb_tableS5(7) zeros(1,13)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,7) min_data_mean_HmAb_tableS5(8) zeros(1,12)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,8) min_data_mean_HmAb_tableS5(9) zeros(1,11)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,9) min_data_mean_HmAb_tableS5(10) zeros(1,10)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,10) min_data_mean_HmAb_tableS5(11) zeros(1,9)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,11) min_data_mean_HmAb_tableS5(12) zeros(1,8)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,12) min_data_mean_HmAb_tableS5(13) zeros(1,7)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,13) min_data_mean_HmAb_tableS5(14) zeros(1,6)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,14) min_data_mean_HmAb_tableS5(15) zeros(1,5)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,15) min_data_mean_HmAb_tableS5(16) zeros(1,4)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,16) min_data_mean_HmAb_tableS5(17) zeros(1,3)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,17) min_data_mean_HmAb_tableS5(18) zeros(1,2)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,18) min_data_mean_HmAb_tableS5(19) zeros(1,1)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);
bar(1:20, [zeros(1,19) min_data_mean_HmAb_tableS5(20) zeros(1,0)], 0.4, 'FaceColor',darkgray(1,:),'FaceAlpha',0.7);

set(gca,'XTick',1:20,'XTickLabel',...
    {'AR3A','AR3B','AR3C',...
    'AR3D','HCV1','HC33-32','HC-1','1:7','A8','e137',...
    'CBH-5','HC84-1','HC84-21', 'HC84-22', 'HC84-23',...
    'HC84-25','HC84-27','HC33-8','HC33-29','CBH-2'});
xtickangle(45)
ylabel('Minimum escape time')
xlim([0.5 20.5])
post_proc_fig


%Heatmap

thresh_range = 80:20:140;
for mm = 1:length(thresh_range)
    for kk = 1:length(data_mean_HmAb_tableS5)
        no_of_sites_less_100_HmAb_tableS5(kk,mm) = sum(data_mean_HmAb_tableS5{kk}<thresh_range(mm));
    end
end

labels_cell{2} = {'80','100','120','140'};

labels_cell{1} = {'AR3A','AR3B','AR3C',...
    'AR3D','HCV1','HC33-32','HC-1','1:7','A8','e137',...
    'CBH-5','HC84-1','HC84-21', 'HC84-22', 'HC84-23',...
    'HC84-25','HC84-27','HC33-8','HC33-29','CBH-2'};

% figure_heatmap(no_of_sites_less_100_HmAb40.','BuPu','',{'','Escape time'},[0 8],labels_cell)

subplot(2,1,2)
h = heatmap(labels_cell{1},labels_cell{2},no_of_sites_less_100_HmAb_tableS5.');
set(gca,'Colormap',color_scheme_BuPu,'ColorLimits',[0 14],'GridVisible','off','FontSize',10,'ColorbarVisible','off');
% h = set(gca,'Colormap',color_scheme,'ColorLimits',limits_data);
% xlabel(text_labels{1})
ylabel('Escape time threshold, tau')

