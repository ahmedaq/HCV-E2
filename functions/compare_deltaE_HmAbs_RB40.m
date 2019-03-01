function compare_deltaE_HmAbs_RB40(mean_escape_time)

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

%% Pierce2016, Table S1 [red i.e., <=40%]

HmAb_Pierce2016 = [];

%CBH-4D
HmAb_Pierce2016{1} = [452 459 494 495 496 497 499 502 504 506:509 511 514 516 537:539 540 542:545 547:552 554:556 558 559 561 562 564 565 569 583 584 585 592:594 597 598 600 602 603 604 607:611 612 614 617:621 623 624:626 627:633 634 636:638 640 642:644];
%CBH-4G
HmAb_Pierce2016{2} = [452 459 494 495 496 497 499 503 505:510 511 513 516 517 537 538 539 540 542:545 547:552 554:556 558 559 561 562 564 565 569 576 578 583 584 587 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624:627 629 631:633 638 640 642:644];
%CBH-4B
HmAb_Pierce2016{3} = [494 497 503 505:511 517 537:539 540 542:545 547 549:552 554 555 556 559 561 562 564 565 569 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624:627 629 630:633 637 638 640 642:644];
%CBH-20
HmAb_Pierce2016{4} = [452 459 494 497 503 505:511 516 517 537:539 540 542:545 547:552 554:556 558 559 561 562 564 565 569 578 581 583:585 592:594 597 598 600 602 603 604 607 608:611 614 617:620 621 623 624 627 629 631 632 633 638 640 642:644];
%CBH-21
HmAb_Pierce2016{5} = [452 459 494 497 503 505:511 517 537:539 540 542:545 547:552 554 555 556 559 561 562 564 565 569 581 583:585 592:594 597 598 600 602 603 604 607:611 614 617:621 623 624 625 627 629 631:633 638 640 642:644];
%CBH-22
HmAb_Pierce2016{6} = [452 459 494 497 503 505:511 517 537:539 540 542:545 547:552 554:556 559 561 562 564 565 566 569 578 581 583:585 592:594 597 598 600 602 603 604 607:611 614 617:621 623 624 625 627 629 631:633 638 640 642:644];

%HC-1
HmAb_Pierce2016{7} = [426 428:430 494 497 503:509 529 530 535 537 539 552 554 559 564 602 603 607 608 610 611 614 617:619 621 624 640 642:644];
%HC-11
HmAb_Pierce2016{8} = [425:429 436:438 442 443 446 459 494 496 497 502:509 511 516 518 520 523 529 530 535 537 539 540 543 547 550:552 554:556 558 559 562 564 565 569 594 597 600 602:604 607:611 614 617:621 623:626 638 640 642:644];
%A27
HmAb_Pierce2016{9} = [424:429 432 437 438 441 442 459 494 496 497 499 502:509 511 516 518:520 523 529 530 531 535:540 550:552 554:556 558 559 561 562 564 565 600 602 603 607:611 614 616:619 621 623 624 625 638 640 642:644];

%CBH-23
HmAb_Pierce2016{10} = [494 497 505:509 537 539 549 552 554 564 607 611 614 617 618 621 643 644]; 
%CBH-7
HmAb_Pierce2016{11} = [494 497 504 506:509 537 539 544 547 549 550 552 554 559 564 602 603 607 611 614 617 618 621 624 626 640 642:644];

%HC84-20
HmAb_Pierce2016{12} = [429 437 441 442 443 459 494 496 497 499 502:509 511 515 516 519 520 523 535 537 539 540 550:552 554:556 558 559 562 564 565 600 602 603 607 608 609 610 611 613 614 616:620 621 623 624 638 640 642 643 644];
%HC84-24
HmAb_Pierce2016{13} = [429 442 443 459 494 496 497 499 502:509 511 516 537 539 540 550:552 554 556 559 564 565 602 603 607 608:611 614 617:619 621 623 624 638 640 642 643 644];
%HC84-26
HmAb_Pierce2016{14} = [429 441 442 446 459 494 497 499 502:509 511 537 539 540 551 552 554 559 564 565 602 603 607 608:611 614 616:619 621 624 638 640 642 643 644];

%HC33-1
HmAb_Pierce2016{15} = [413 418 420 652];
%HC33-4
HmAb_Pierce2016{16} = [408 413 420 652];
% %CD81_bs
% HmAb_Pierce2016{17} = [420 421 424 427 430 432 436:438 440:443 523 526 527 529 530 535 540 549 550 613 614 616:618];


for kk = 1:length(HmAb_Pierce2016)
    data_mean_HmAb20{kk} = mean_escape_time(HmAb_Pierce2016{kk}-383); 
    min_data_mean_HmAb20(kk) = min(data_mean_HmAb20{kk});
end

ylim_min = 0;
ylim_max = 300;

figure
subplot(2,1,1)

% xbars = [6.5 7.5];
% patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
hold on
xbars = [9.5 10.5];
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

%% Pierce2016, Table S1 [red and green i.e., <=40%] v2
% 
% HmAb_Pierce2016 = [];
% 
% %CBH-4D
% HmAb_Pierce2016{1} = [452 459 494 495 496 497 499 502 504 506:509 511 514 516 537:539 540 542:545 547:552 554:556 558 559 561 562 564 565 569 583 584 585 592:594 597 598 600 602 603 604 607:611 612 614 617:621 623 624:626 627:633 634 636:638 640 642:644];
% %CBH-4G
% HmAb_Pierce2016{2} = [452 459 494 495 496 497 499 503 505:510 511 513 516 517 537 538 539 540 542:545 547:552 554:556 558 559 561 562 564 565 569 576 578 583 584 587 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624:627 629 631:633 638 640 642:644];
% %CBH-4B
% HmAb_Pierce2016{3} = [494 497 503 505:511 517 537:539 540 542:545 547 549:552 554 555 556 559 561 562 564 565 569 592:594 597 598 600 602 603 607:611 614 617:619 621 623 624:627 629 630:633 637 638 640 642:644];
% %CBH-20
% HmAb_Pierce2016{4} = [452 459 494 497 503 505:511 516 517 537:539 540 542:545 547:552 554:556 558 559 561 562 564 565 569 578 581 583:585 592:594 597 598 600 602 603 604 607 608:611 614 617:620 621 623 624 627 629 631 632 633 638 640 642:644];
% %CBH-21
% HmAb_Pierce2016{5} = [452 459 494 497 503 505:511 517 537:539 540 542:545 547:552 554 555 556 559 561 562 564 565 569 581 583:585 592:594 597 598 600 602 603 604 607:611 614 617:621 623 624 625 627 629 631:633 638 640 642:644];
% %CBH-22
% HmAb_Pierce2016{6} = [452 459 494 497 503 505:511 517 537:539 540 542:545 547:552 554:556 559 561 562 564 565 566 569 578 581 583:585 592:594 597 598 600 602 603 604 607:611 614 617:621 623 624 625 627 629 631:633 638 640 642:644];
% 
% %HC-1
% HmAb_Pierce2016{7} = [426 428:430 494 497 503:509 529 530 535 537 539 552 554 559 564 602 603 607 608 610 611 614 617:619 621 624 640 642:644];
% %HC-11
% HmAb_Pierce2016{8} = [425:429 436:438 442 443 446 459 494 496 497 502:509 511 516 518 520 523 529 530 535 537 539 540 543 547 550:552 554:556 558 559 562 564 565 569 594 597 600 602:604 607:611 614 617:621 623:626 638 640 642:644];
% %A27
% HmAb_Pierce2016{9} = [424:429 432 437 438 441 442 459 494 496 497 499 502:509 511 516 518:520 523 529 530 531 535:540 550:552 554:556 558 559 561 562 564 565 600 602 603 607:611 614 616:619 621 623 624 625 638 640 642:644];
% 
% %CBH-23
% HmAb_Pierce2016{10} = [494 497 505:509 537 539 549 552 554 564 607 611 614 617 618 621 643 644]; 
% %CBH-7
% HmAb_Pierce2016{11} = [494 497 504 506:509 537 539 544 547 549 550 552 554 559 564 602 603 607 611 614 617 618 621 624 626 640 642:644];
% 
% %HC84-20
% HmAb_Pierce2016{12} = [429 437 441 442 443 459 494 496 497 499 502:509 511 515 516 519 520 523 535 537 539 540 550:552 554:556 558 559 562 564 565 600 602 603 607 608 609 610 611 613 614 616:620 621 623 624 638 640 642 643 644];
% %HC84-24
% HmAb_Pierce2016{13} = [429 442 443 459 494 496 497 499 502:509 511 516 537 539 540 550:552 554 556 559 564 565 602 603 607 608:611 614 617:619 621 623 624 638 640 642 643 644];
% %HC84-26
% HmAb_Pierce2016{14} = [429 441 442 446 459 494 497 499 502:509 511 537 539 540 551 552 554 559 564 565 602 603 607 608:611 614 616:619 621 624 638 640 642 643 644];
% 
% %HC33-1
% HmAb_Pierce2016{15} = [413 418 420 652];
% %HC33-4
% HmAb_Pierce2016{16} = [408 413 420 652];
% % %CD81_bs
% % HmAb_Pierce2016{17} = [420 421 424 427 430 432 436:438 440:443 523 526 527 529 530 535 540 549 550 613 614 616:618];
% 
% binding_residues_all_40 = unique([HmAb_Pierce2016{:}]);
% 
% G = [];
% data = [];
% 
% for kk = 1:length(HmAb_Pierce2016)
%     G = [G kk*ones(1,length(HmAb_Pierce2016{kk}))];
%     data = [data dE2(HmAb_Pierce2016{kk}-383)];
% end
% 
% set(0,'DefaultAxesFontName','Arial')
% set(0,'DefaultTextFontName','Arial')
% set(0,'DefaultAxesFontSize',10)
% set(0,'DefaultTextFontSize',10)
% 
% box_lineWidth = 0.5;
% box_widths_value = 0.3;
% % box_color = [repmat(blue,6,1);red;repmat(blue,2,1);repmat(red,3,1);repmat(blue,2,1);red;repmat(blue,2,1)];
% % box_color = [repmat(green,6,1);repmat(purple,3,1);repmat(orange,2,1);repmat(yellow,3,1);repmat(brown,2,1)];
% box_color = [repmat(color_scheme_npg(1,:),6,1);repmat(color_scheme_npg(2,:),3,1);...
%     repmat(color_scheme_npg(3,:),2,1);repmat(color_scheme_npg(4,:),3,1);repmat(color_scheme_npg(5,:),2,1)];
% box_color_transparency = 0.7; %faceAlpha
% median_lineWidth = 1;
% median_color = 'k';
% whisker_value = 1.5;
% outlier_marker = 'o';
% outlier_markerSize = 7;
% outlier_marker_edgeWidth = 0.1;
% outlier_marker_edgeColor = 'k';
% outlier_jitter_value = 0.75;
% label_xaxis_data = {'CBH-4D','CBH-4G','CBH-4B','CBH-20','CBH-21','CBH-22',...
%     'HC-1','HC-11','A27','CBH-23','CBH-7','HC84-20','HC84-24','HC84-26',...
%     'HC33-1','HC33-4'};
% text_ylabel = 'Fitness cost, \DeltaE_i';
% text_xlabel = '';%'HmAb';
% text_title = '';%Pierce2016, Table S1, Binding <=20%';
% label_orientation_choice = 'horizontal'; %'inline'
% ylim_min = -1;
% ylim_max = 13;
% savefig = 0;
% savefig_name = 'comparison_dE_E2_abs_Pierce2016_TableS1_lt40_v2';
% fig_width_cm = 10;%16;
% fig_height_cm = 4;%6
% 
% figure
% xbars = [9.5 10.5];
% patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% xbars = [14.5 15.5];
% patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% 
% figure_boxplot(data,G,...
%     box_lineWidth,box_widths_value,box_color,box_color_transparency,...
%     median_lineWidth,median_color,...
%     whisker_value,...
%     outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
%     label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
%     ylim_min,ylim_max,...
%     savefig,savefig_name,fig_width_cm,fig_height_cm);
% 
% %% Pierce2016, Table S5, without non-neutralizing AR1A and AR1B v2 (without those in Table S1)
% %% Also including other HmAbs 
% % [Owsianka2008]
% thresh_binding_imp = 25; %<%
% [dataInput1] = xlsread('Binding_CHB5_Owsianka2008.xlsx');
% CBH5 =  dataInput1(dataInput1(:,2)<thresh_binding_imp,1);  
% CBH5 =  unique([CBH5.' 494 497 614 617:619 621 624 502:505 507:509]); %extra sites in [lacob2008]
% [dataInput2] = xlsread('Binding_CHB7_Owsianka2008.xlsx');
% CBH7 =  dataInput2(dataInput2(:,2)<thresh_binding_imp,1);
% CBH7 = unique([CBH7.' 494 497 614 617:619 621 624]); %extra sites in [lacob2008]
% 
% HmAb_Pierce2016_tableS5 = [];
% 
% HmAb_Pierce2016_tableS5{1} = [424 523 525 530 535 538 540];%AR3A
% HmAb_Pierce2016_tableS5{2} = [412, 416, 418, 423, 424, 523, 525, 530, 535, 540]; %[424 530 535]; %AR3B
% HmAb_Pierce2016_tableS5{3} = [424, 488, 523, 525, 530, 535, 538, 540];%[424 525 530 535 540]; %AR3C
% HmAb_Pierce2016_tableS5{4} = [412, 424, 523, 530, 535];%[424 530]; %AR3D
% HmAb_Pierce2016_tableS5{5} = [412:423]; %HCV1
% % HmAb_Pierce2016_tableS5{6} = [413 418 420]; %HC33.1
% HmAb_Pierce2016_tableS5{6} = [413 418 420]; %HC33.32
% HmAb_Pierce2016_tableS5{7} = [525 530 535]; %HC-1
% HmAb_Pierce2016_tableS5{8} = CBH5; %[Owsianka2008,lacob2008]
% % HmAb_Pierce2016_tableS5{10} = CBH7; %[Owsianka2008,lacob2008]
% HmAb_Pierce2016_tableS5{9} = [523 526 527 529 530 535]; %1:7
% HmAb_Pierce2016_tableS5{10} = [523 526 527 529 530 535]; %A8
% HmAb_Pierce2016_tableS5{11} = [416 420 529 530 535]; %e137 [Sabo2011]
% 
% HmAb_Pierce2016_tableS5{12} = [441 442]; %HC84.1
% % HmAb_Pierce2016_tableS5{15} = [441:443 613 616]; %HC84.20
% HmAb_Pierce2016_tableS5{13} = [441:443]; %HC84.21
% HmAb_Pierce2016_tableS5{14} = [420 428 437 441:443 616]; %HC84.22
% HmAb_Pierce2016_tableS5{15} = [420 428 437 441:443 616]; %HC84.23
% % HmAb_Pierce2016_tableS5{19} = [442:443]; %HC84.24
% HmAb_Pierce2016_tableS5{16} = [441 442 616]; %HC84.25
% % HmAb_Pierce2016_tableS5{21} = [441 442 616]; %HC84.26
% HmAb_Pierce2016_tableS5{17} = [441:443 446 616]; %HC84.27
% % HmAb_Pierce2016_tableS5{23} = [408 413 420]; %HC33.4
% HmAb_Pierce2016_tableS5{18} = [408 413 418 420]; %HC33.8
% HmAb_Pierce2016_tableS5{19} = [408 413 418 420]; %HC33.29
% % HmAb_Pierce2016_tableS5{26} = [425 428 436:438 442 443 530 535]; %HC-11
% HmAb_Pierce2016_tableS5{20} = [437 439 530 535 523 431]; %added residues from [Keck2008] %[437 439 530 535]; %CBH-2
% 
% %
% 
% G = [];
% data = [];
% 
% for kk = 1:length(HmAb_Pierce2016_tableS5)
%     G = [G kk*ones(1,length(HmAb_Pierce2016_tableS5{kk}))];
%     data = [data dE2(HmAb_Pierce2016_tableS5{kk}-383)];
% end
% 
% 
% set(0,'DefaultAxesFontName','Arial')
% set(0,'DefaultTextFontName','Arial')
% set(0,'DefaultAxesFontSize',10)
% set(0,'DefaultTextFontSize',10)
% 
% box_lineWidth = 0.5;
% box_widths_value = 0.3;
% box_color = darkgray;%[repmat(red,11,1);repmat(blue,9,1)];
% box_color_transparency = 0.5; %faceAlpha
% median_lineWidth = 1;
% median_color = 'k';
% whisker_value = 1.5;
% outlier_marker = 'o';
% outlier_markerSize = 7;
% outlier_marker_edgeWidth = 0.1;
% outlier_marker_edgeColor = 'k';
% outlier_jitter_value = 0.75;
% label_xaxis_data = {'AR3A','AR3B','AR3C',...
%     'AR3D','HCV1','HC33-32','HC-1','CBH-5','1:7','A8','e137',...
%     'HC84-1','HC84-21', 'HC84-22', 'HC84-23',...
%     'HC84-25','HC84-27','HC33-8','HC33-29','CBH-2'};
% text_ylabel = 'Fitness cost, \DeltaE_i';
% text_xlabel = '';%'HmAb';
% text_title = '';%Other HmAbs [Pierce2016,Table S5],[Keck2008],[Owsianka2008],[Iacob2008],[Perotti2008]';
% label_orientation_choice = 'horizontal'; %'inline'
% ylim_min = -1;
% ylim_max = 13;
% savefig = 0;
% savefig_name = 'comparison_dE_E2_abs_Pierce2016_TableS5_v2';
% fig_width_cm = 10;%16;
% fig_height_cm = 4;%6;
% 
% figure
% xbars = [0.525 11.5];
% patch([xbars(1) xbars(1), xbars(2) xbars(2)], [ylim_min+.05 ylim_max ylim_max ylim_min+.05], [0.85 0.85 0.85], ...
%     'EdgeColor','w','LineWidth',0.1)
% hold on
% 
% figure_boxplot(data,G,...
%     box_lineWidth,box_widths_value,box_color,box_color_transparency,...
%     median_lineWidth,median_color,...
%     whisker_value,...
%     outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
%     label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
%     ylim_min,ylim_max,...
%     savefig,savefig_name,fig_width_cm,fig_height_cm);