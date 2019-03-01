function compare_deltaE_HmAbs(mean_escape_time)

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

%% Pierce2016, Table S1 [red i.e., <=20%] v2

HmAb_Pierce2016 = [];

%CBH-4D
HmAb_Pierce2016{1} = [494 497 502 504 506:509 511 537 539 540 542:545 547 549:552 554 556 559 561 562 564 565 584 585 592 594 598 600 602 603 607:611 614 617:619 621 623 624 626 627 629:633 638 640 642:644];
%CBH-4G
HmAb_Pierce2016{2} = [494 497 503 505:509 511 517 537 539 540 542:545 547 549:552 554 556 559 561 564 565 584 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 629 631:633 638 640 642:644];
%CBH-4B
HmAb_Pierce2016{3} = [494 497 503 505:509 537 539 540 550:552 554 559 561 564 565 600 602 603 607:611 614 617:619 621 624 627 629 631:633 638 640 642:644];
%CBH-20
% HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 611 614 617:619 621 623 624 627 631 632 638 640 642:644];
%included the residues with close to threshold (RB = 21,22) as well
HmAb_Pierce2016{4} = [494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592 594 597 598 600 602 603 607 608 610 611 614 617:619 621 623 624 627 631 632 638 640 642:644]; 
%CBH-21
HmAb_Pierce2016{5} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];
%CBH-22
HmAb_Pierce2016{6} = [452 494 497 503 505:509 537 539 540 542:545 547 549:552 554 559 561 564 565 569 581 583:585 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624 627 631:633 638 640 642:644];

%HC-1
HmAb_Pierce2016{7} = [429 494 503:506 508 509 529 530 535 537 539 552 554 559 564 607 611 614 617 644];
%HC-11
HmAb_Pierce2016{8} = [425 428 429 436:438 442 443 494 497 502:504 506:509 511 520 530 535 537 539 550:552 554 556 559 564 565 602 603 607 608 611 614 617:621 624 640 643 644];
%A27
HmAb_Pierce2016{9} = [424:429 437 438 494 497 499 502:504 506:509 511 520 529 530 535 537 539 540 550:552 554 556 559 564 565 602 603 607:611 614 616:619 621 623 624 638 640 642:644];

%CBH-23
HmAb_Pierce2016{10} = [494 508 509 537 539 549 552 554 564 611 614 644]; %same as CBH7
%CBH-7
HmAb_Pierce2016{11} = [494 506 508 509 537 539 549 552 554 564 611 614 617 621 644];

%HC84-20
HmAb_Pierce2016{12} = [429 441 494 497 502:509 511 537 539 552 554 559 564 607 608 611 613 614 616:619 621 640 643 644];
%HC84-24
HmAb_Pierce2016{13} = [429 442 443 494 497 502:509 511 537 539 552 554 559 564 607 608 611 614 617:619 621 643 644];
%HC84-26
HmAb_Pierce2016{14} = [441 442 494 497 502:509 511 537 539 552 554 559 564 603 607 608 611 614 617:619 621 640 643 644];

%HC33-1
HmAb_Pierce2016{15} = [413 418 420];
%HC33-4
HmAb_Pierce2016{16} = [408 413 420];
% %CD81_bs
% HmAb_Pierce2016{17} = [420 421 424 427 430 436:438 440:443 523 527 529 530 535 540 613 614 616:618];


for kk = 1:length(HmAb_Pierce2016)
    data_mean_HmAb20{kk} = mean_escape_time(HmAb_Pierce2016{kk}-383); 
    min_data_mean_HmAb20(kk) = min(data_mean_HmAb20{kk});
end

ylim_min = 0;
ylim_max = 300;

figure
subplot(2,1,1)

xbars = [6.5 7.5];
patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
    'EdgeColor','w','LineWidth',0.1)
hold on
xbars = [9.5 11.5];
patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
    'EdgeColor','w','LineWidth',0.1)
xbars = [14.5 15.5];
patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
    'EdgeColor','w','LineWidth',0.1)

bar(1:16, [zeros(1,0) min_data_mean_HmAb20(1) zeros(1,15)], 0.4, 'FaceColor',color_scheme_npg(1,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,1) min_data_mean_HmAb20(2) zeros(1,14)], 0.4, 'FaceColor',color_scheme_npg(1,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,2) min_data_mean_HmAb20(3) zeros(1,13)], 0.4, 'FaceColor',color_scheme_npg(1,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,3) min_data_mean_HmAb20(4) zeros(1,12)], 0.4, 'FaceColor',color_scheme_npg(1,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,4) min_data_mean_HmAb20(5) zeros(1,11)], 0.4, 'FaceColor',color_scheme_npg(1,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,5) min_data_mean_HmAb20(6) zeros(1,10)], 0.4, 'FaceColor',color_scheme_npg(1,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,6) min_data_mean_HmAb20(7) zeros(1,9)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,7) min_data_mean_HmAb20(8) zeros(1,8)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,8) min_data_mean_HmAb20(9) zeros(1,7)], 0.4, 'FaceColor',color_scheme_npg(2,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,9) min_data_mean_HmAb20(10) zeros(1,6)], 0.4, 'FaceColor',color_scheme_npg(3,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,10) min_data_mean_HmAb20(11) zeros(1,5)], 0.4, 'FaceColor',color_scheme_npg(3,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,11) min_data_mean_HmAb20(12) zeros(1,4)], 0.4, 'FaceColor',color_scheme_npg(4,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,12) min_data_mean_HmAb20(13) zeros(1,3)], 0.4, 'FaceColor',color_scheme_npg(4,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,13) min_data_mean_HmAb20(14) zeros(1,2)], 0.4, 'FaceColor',color_scheme_npg(4,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,14) min_data_mean_HmAb20(15) zeros(1,1)], 0.4, 'FaceColor',color_scheme_npg(5,:),'FaceAlpha',0.7);
bar(1:16, [zeros(1,15) min_data_mean_HmAb20(16) zeros(1,0)], 0.4, 'FaceColor',color_scheme_npg(5,:),'FaceAlpha',0.7);

set(gca,'XTick',1:16,'XTickLabel',...
    {'CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
    'HC-1','HC-11','A27','CBH-23','CBH-7','HC84-20','HC84-24','HC84-26',...
    'HC33-1','HC33-4'});
xtickangle(45)
ylabel('Minimum escape time')
xlim([0.5 16.5])
post_proc_fig


%Heatmap

thresh_range = 80:20:140;
for mm = 1:length(thresh_range)
    for kk = 1:length(data_mean_HmAb20)
        no_of_sites_less_100_HmAb20(kk,mm) = sum(data_mean_HmAb20{kk}<=thresh_range(mm));
    end
end

labels_cell{2} = {'80','100','120','140'};
labels_cell{1} = {'CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
    'HC-1','HC-11','A27','CBH-23','CBH-7','HC84-20','HC84-24','HC84-26',...
    'HC33-1','HC33-4'};

% figure_heatmap(no_of_sites_less_100_HmAb20.','BuPu','',{'','Escape time'},[0 8],labels_cell)

subplot(2,1,2)
heatmap(labels_cell{1},labels_cell{2},no_of_sites_less_100_HmAb20.');
set(gca,'Colormap',color_scheme_BuPu,'ColorLimits',[0 10],'GridVisible','off','FontSize',10,'ColorbarVisible','off');
ylabel('Escape time threshold, tau')
% colorbar('Location','north')
