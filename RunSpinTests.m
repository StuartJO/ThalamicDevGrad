PC1 = load('VERT_WEI_PC1_Thr100_Results.mat');
PC2 = load('VERT_WEI_PC2_Thr100_Results.mat');
PC3 = load('VERT_WEI_PC3_Thr100_Results.mat');

load('VERTWEI_Thr100.mat');
%medwallmask = PC1.medwallmask;    

SurfData = NaN(length(medwallmask),1);

PREvsTERM_diff = mean(PRETERM_Cort_TC_norm)-mean(TERM_Cort_TC_norm);

[~,~,~,stats] = ttest2(PRETERM_Cort_TC_norm,mean(TERM_Cort_TC_norm));

PREvsTERM_diff = stats.tstat;

PREvsTERM_diff_thal = mean(PRETERM_Thal_TC_norm)-mean(TERM_Thal_TC_norm);

[~,~,~,stats] = ttest2(PRETERM_Thal_TC_norm,mean(TERM_Thal_TC_norm));

PREvsTERM_diff_thal = stats.tstat;

SurfData(medwallmask) = PREvsTERM_diff;



DATA(:,1) = TrainCoeffz(:,1);
DATA(:,2) = TrainCoeffz(:,2);
DATA(:,3) = TrainCoeffz(:,3);

DATA(:,4) = PC1.PC_age_corr;
DATA(:,5) = PC2.PC_age_corr;
DATA(:,6) = PC3.PC_age_corr;

DATA_LABEL = {'PC1 template loading','PC2 template loading','PC3 template loading','PC1 age loading corr.','PC2 age loading corr.','PC3 age loading corr.'};
DATA2_LABEL = {'PC1 template score','PC2 template score','PC3 template score','PC1 age score corr.','PC2 age score corr.','PC3 age score corr.'};

DATA2(:,1) = TrainScorez(:,1);
DATA2(:,2) = TrainScorez(:,2);
DATA2(:,3) = TrainScorez(:,3);

DATA2(:,4) = PC1.PC_thal_age_corr;
DATA2(:,5) = PC2.PC_thal_age_corr;
DATA2(:,6) = PC3.PC_thal_age_corr;

for i = 1:6

InDATA = NaN(length(medwallmask),1);
InDATA(medwallmask) = DATA(:,i);

p_perm = perm_sphere_p(InDATA,SurfData,spins','Pearson');

%figure('Position', [291 391 560 420])
subplot(2,3,i)
        
if p_perm == 0
    pval_format = '{\itp_{spin}} < .001 ';   
else
    pval_format = ['{\itp_{spin}} = ', num2str(round(p_perm, 4))];
end
s = scatterfit(InDATA, SurfData, 40,lines(1),pval_format,2);
%title(['{\itr} = ', num2str(neuromap_corrs.corr(SigInd), 3), ', ', pval_format, ''], 'FontWeight', 'normal');
% xlabel(xlabel_var)
% ylabel(neuromap_corrs.name{SigInd})
% colormap(turbo)
xlabel(DATA_LABEL{i})
ylabel('Cortical preterm-term diff.')

set(gca, 'FontSize', 16)
s.MarkerEdgeColor = [0 0 0];

end

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_cort_diff_scatts.png'],'-dpng','-r300')

figure 
for i = 1:6
    InDATA =  DATA2(:,i);
    subplot(2,3,i)
    s = scatterfit(InDATA,PREvsTERM_diff_thal,40,lines(1),[],2);
    xlabel(DATA2_LABEL{i})
ylabel('Thalamic preterm-term diff.')
    set(gca, 'FontSize', 16)
s.MarkerEdgeColor = [0 0 0];
end

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_thal_diff_scatts.png'],'-dpng','-r300')
