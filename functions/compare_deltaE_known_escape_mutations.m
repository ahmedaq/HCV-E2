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
    [391 394 401 415 417 434 444 608];

%Escape mutation [Keck2008]
polymorphisms_associated_with_neutralization_resistance4 = 431;

%Escape mutation from HC33-4 [Keck2016]
polymorphisms_associated_with_neutralization_resistance5 = 408;

%Escape mutations in HVR1 [Kato1993]
polymorphisms_associated_with_neutralization_resistance6 = ...
    [384 386 388 390 393 395 399 401 391 394 396:400 410 407 402:405];


%%

polymorphisms_associated_with_neutralization_resistance = unique([...
    polymorphisms_associated_with_neutralization_resistance1 ...
    polymorphisms_associated_with_neutralization_resistance2 ...
    polymorphisms_associated_with_neutralization_resistance3 ...
    polymorphisms_associated_with_neutralization_resistance4 ...
    polymorphisms_associated_with_neutralization_resistance5 ...
    polymorphisms_associated_with_neutralization_resistance6]);

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
box_color = [blue;red];
box_color_transparency = 0.7; %faceAlpha
median_lineWidth = 1;
median_color = 'k';
whisker_value = 1.5;
outlier_marker = 'o';
outlier_markerSize = 4;
outlier_marker_edgeWidth = 0.1;
outlier_marker_edgeColor = 'w';
outlier_jitter_value = 0;
label_xaxis_data = {'Escape','Remaining'};
text_ylabel = 'Escape time (generations)';
text_xlabel = '';
text_title = '';%'E2-escape mutations [Keck2009],[Morin2012],[Bailey2015]';
label_orientation_choice = 'horizontal'; %'horizontal'
ylim_min = 0;
ylim_max = 1000;
savefig = 0;
savefig_name = 'escape_mutations';
fig_width_cm = 4;
fig_height_cm = 5;

figure;
figure_boxplot(data,G,...
    box_lineWidth,box_widths_value,box_color,box_color_transparency,...
    median_lineWidth,median_color,...
    whisker_value,...
    outlier_marker,outlier_markerSize,outlier_marker_edgeWidth,outlier_marker_edgeColor,outlier_jitter_value,...
    label_xaxis_data,text_ylabel,text_xlabel,text_title,label_orientation_choice,...
    ylim_min,ylim_max,...
    savefig,savefig_name,fig_width_cm,fig_height_cm);


P = ranksum(dE2(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),dE2(polymorphisms_associated_with_neutralization_resistance-383),'tail','right');
fprintf('\nP = %.1e, Mann-Whitney test\n',P)
% [h,p] = kstest2(dE2(setdiff(1:L,polymorphisms_associated_with_neutralization_resistance-383)),dE2(polymorphisms_associated_with_neutralization_resistance-383),'tail','smaller')

%% AUC curve for a classifier based on mean escape time

% set(0,'DefaultTextFontSize',8)
% set(0,'DefaultAxesFontSize',8)
ls = 363;

correct_labels_exp = zeros(ls,1);
correct_labels_exp(polymorphisms_associated_with_neutralization_resistance-383) = 1;

[X,Y,~,AUC] = perfcurve(correct_labels_exp,dE2(:),0);

figure;
plot(X,Y,'LineWidth',1.5,'Color',blue)
hold on
% title(sprintf('AUC = %.3f',AUC))
set(gca, ...
    'Box'         , 'on'     , ...
    'TickDir'     , 'in'     )
xlabel('False positive rate'); ylabel('True positive rate')
text(0.6, 0.2, sprintf('AUC = %.2f',AUC),'FontSize',8,'FontWeight','bold')

% set(gcf,'PaperUnits','centimeters','PaperPosition',[0 0 8 6])
% print('AUC_escape_mutants','-dpng','-r600')


% [X,Y] = perfcurve(correct_labels_exp,mean_escape_time(:),0,'nboot',1000,'xvals','all');
% errorbar(X,Y(:,1),Y(:,1)-Y(:,2),Y(:,3)-Y(:,1),'Color',blue); % plot errors

%% designing threshold based on max accuracy/F1 score

threshold = 0:1:round(max(dE2));

Fscore = zeros(1,length(threshold));
accuracy = zeros(1,length(threshold));
specificity = zeros(1,length(threshold));
tpr = zeros(1,length(threshold));
fpr = zeros(1,length(threshold));
MCC = zeros(1,length(threshold));

class = 1;

%%% Using classification_performance.m
for kk = 1:length(threshold)
    
    thresh = threshold(kk);
    
    group_vec = correct_labels_exp;
    group_vec_hat = double(dE2(:)<=thresh);
    
    
    
    [Fscore_out, accuracy_out, tpr_out, specificity_out, ~, ~, TP_out, TN_out, FP_out, FN_out] = ...
        classification_performance(group_vec,group_vec_hat);
    
    
    Fscore(kk) = Fscore_out(class);
    accuracy(kk) = accuracy_out(class);
    tpr(kk) = tpr_out(class);
    specificity(kk) = specificity_out(class);
    
    
    TP = TP_out(class);
    FP = FP_out(class);
    TN = TN_out(class);
    FN = FN_out(class);
    
    MCC(kk) = ((TP*TN)- (FP*FN)) / (sqrt( (TP+FP)*(TP+FN)*(TN+FP)*(TN+FN) ));
    if isnan(MCC(kk))
        MCC(kk) = (TP*TN)- (FP*FN);
    end
    
end
fpr = 1-specificity;


[max_Fscore,indx_max_Fscore] = max(Fscore)
thresh_max_Fscore = threshold(indx_max_Fscore)

[max_acc,indx_max_acc] = max(accuracy)
thresh_max_accuracy = threshold(indx_max_acc)

% [max_tpr,indx_max_tpr] = max(tpr)
% thresh_max_tpr = threshold(indx_max_tpr)
% 
% [min_fpr,indx_min_fpr] = min(fpr)
% thresh_min_fpr = threshold(indx_min_fpr)

[max_MCC,indx_max_MCC] = max(MCC)
thresh_max_MCC = threshold(indx_max_MCC)


%%
% set(0,'DefaultTextFontSize',8)
% set(0,'DefaultAxesFontSize',8)

figure;
hold on
% plot(threshold,accuracy,'Color',red,'LineWidth',1.5)
plot(threshold,Fscore,'Color',orange,'LineWidth',1.5)
plot(threshold,MCC,'Color',brown,'LineWidth',1.5)
% h = legend('Accuracy','F1 score','MCC'); legend boxoff
h = legend('F1 score','MCC'); legend boxoff
xlabel('Escape time score cut-off (generations)')
ylabel('Classification metric score')

plot(thresh_max_Fscore*ones(1,11),0:0.1:1,'--','Color','k','LineWidth',1)
% plot(thresh_max_MCC*ones(1,11),0:0.1:1,':','Color','k','LineWidth',.5)

post_proc_fig
set(gca,'box','on','TickDir','in')
 h.String = {'F1 score'  'MCC'};
text(thresh_max_Fscore-30,1.04,sprintf('%d',thresh_max_Fscore),'Color','k','FontWeight','bold')

