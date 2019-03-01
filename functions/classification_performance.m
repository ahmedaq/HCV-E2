function [Fscore, accuracy, sensitivity, specificity, precision, recall, TP, TN, FP, FN] = ...
    classification_performance(group_vec,group_vec_hat)


stats_clustering = confusionmatStats(group_vec,group_vec_hat);
Fscore = (stats_clustering.Fscore);
accuracy = (stats_clustering.accuracy);
sensitivity = (stats_clustering.sensitivity);
specificity = (stats_clustering.specificity);
precision = (stats_clustering.precision);
recall = (stats_clustering.recall);
TP = stats_clustering.TP;
TN = stats_clustering.TN;
FP = stats_clustering.FP;
FN = stats_clustering.FN;