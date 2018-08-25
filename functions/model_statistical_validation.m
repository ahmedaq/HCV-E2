function [freqCountNumMutPerSeq_actualMSA,freqCountNumMutPerSeq_MCMC] = model_statistical_validation(msa_aa_ex,weight,samples_MCMC_double,phi_curr,protein_length_aa)
  
% Code for statistical validation of the model
% 
% Written by: Ahmed Abdul Quadeer 
% Last updated: 2018-04-07

%%

cross_prod_site = msa_aa_ex'*diag(weight)*msa_aa_ex/sum(weight);

double_mutant_sample = (samples_MCMC_double')*samples_MCMC_double/size(samples_MCMC_double,1);
single_mutant_sample = mean(samples_MCMC_double);

cross_nondiag = cross_prod_site - diag(diag(cross_prod_site));
total_nondiag = double_mutant_sample - diag(diag(double_mutant_sample));

cross_diag = diag(cross_prod_site);
total_diag = diag(double_mutant_sample);

%%%%%%%%%

set(0,'DefaultTextFontSize',10)
set(0,'DefaultAxesFontSize',10)

color = [0.6000    0.6000    0.6000];%gray
markersize = 3;

figure
plot(cross_diag(:),total_diag(:),'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize);hold on
max_value = max(max(cross_diag(:),total_diag(:)));
plot(0:max_value/20:max_value,0:max_value/20:max_value,'k')
xlabel('Single mutant probability (MSA)')
ylabel('Single mutant probability')
% title(['\epsilon_p = ' num2str(sum_error_diag)])
axis([0 0.5 0 0.5])
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'XTick'       , 0:0.1:0.5, ...
    'YTick'       , 0:0.1:0.5, ...
    'LineWidth'   , .5       );

figure
plot(cross_nondiag(:),total_nondiag(:),'o','MarkerEdgeColor',color,'MarkerFaceColor',color,'MarkerSize',markersize);hold on
max_value = max(max(cross_nondiag(:),total_nondiag(:)));
plot(0:max_value/20:max_value,0:max_value/20:max_value,'k')
xlabel('Double mutant probability (MSA)')
ylabel('Double mutant probability')
% title( ['\epsilon_c = ' num2str(sum_error_nondiag) ] )
set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'off'      , ...
    'YMinorTick'  , 'off'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'XTick'       , 0:0.1:0.5, ...
    'YTick'       , 0:0.1:0.5, ...
    'LineWidth'   , .5       );


[xMutMtxCell_actualMSA, allMutMtx_actualMSA] = mutCountMSA(msa_aa_ex, phi_curr, protein_length_aa);
[xMutMtxCell_MCMC, allMutMtx_MCMC] = mutCountMSA(samples_MCMC_double, phi_curr, protein_length_aa);

[freqCountNumMutPerSeq_actualMSA numSeq_actualMSA] = plotXMutFreqPMF(allMutMtx_actualMSA, protein_length_aa);
[freqCountNumMutPerSeq_MCMC numSeq_MCMC] = plotXMutFreqPMF(allMutMtx_MCMC, protein_length_aa);

pdf_actualMSA = freqCountNumMutPerSeq_actualMSA/numSeq_actualMSA./(sum(freqCountNumMutPerSeq_actualMSA/numSeq_actualMSA));
pdf_MCMC = freqCountNumMutPerSeq_MCMC/numSeq_MCMC./(sum(freqCountNumMutPerSeq_MCMC/numSeq_MCMC));


figure
semilogy(0:protein_length_aa, smooth(pdf_actualMSA),'Color','k','LineWidth',1)
hold on
semilogy(0:protein_length_aa, smooth(pdf_MCMC),'Color',color,'LineWidth',1)
% title(sprintf('Distribtion of mutations in %s',heading))
h = legend('MSA','Model');
h.Position = [0.700 0.8524 0.1080 0.0488];
legend boxoff
xlabel('Number of mutations, \itr')
ylabel('Probability of {\itr} mutations')
% no_of_nonzero_mutations = find(pdf_actualMSA ~= 0);
xlim([0 70])
% xlim([-0.1 no_of_nonzero_mutations(end)])
ylim([1e-4 1])

[r,pval] = corr(pdf_MCMC(:),pdf_actualMSA(:),'type','pearson','tail','both');

set(gca, ...
    'Box'         , 'off'     , ...
    'TickDir'     , 'out'     , ...
    'TickLength'  , [.02 .02] , ...
    'XMinorTick'  , 'on'      , ...
    'YMinorTick'  , 'on'      , ...
    'XColor'      , [.1 .1 .1], ...
    'YColor'      , [.1 .1 .1], ...
    'YTick'       , [1e-4 1e-3 1e-2 1e-1 1], ...
    'LineWidth'   , .5       );
