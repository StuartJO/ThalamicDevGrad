
PC1 = load('VERT_WEI_PC1_Thr100_Results.mat');
PC2 = load('VERT_WEI_PC2_Thr100_Results.mat');
PC3 = load('VERT_WEI_PC3_Thr100_Results.mat');


scatterfit

scatter(PC1.PC_thal_age_corr,PC2.PC_thal_age_corr,40,TrainScorez(:,1),'filled')



scatter(TrainScorez(:,1),TrainScorez(:,2),40,PC2.PC_thal_age_corr,'filled')

scatter(TrainScorez(:,1),TrainScorez(:,2),40,PC2.PC_thal_age_corr,'filled')

scatterfit(TrainScorez(:,1),PC1.PC_thal_age_corr,40,'filled')



scatterfit(tiedrank(TrainCoeffz(:,3)),tiedrank(PC3.PC_age_corr),40,'filled');


scatterfit(tiedrank(TrainCoeffz(:,1)),tiedrank(PC1.PC_age_corr),40,'filled');


scatterfit((TrainCoeffz(:,1)),(PC1.PC_age_corr),40,'filled');

for PC = 1:3
figure('Position',[62   300   1300   335])
for i = 1:3
    subplot(1,3,i)
[s,t] = scatterfit(vox_coords(:,i),TrainScorez(:,PC),36);
s.MarkerFaceAlpha=1;
cardDir={'Medial-Lateral','Anterior-posterior','Ventral-Dorsal'};
xlabel(cardDir{i})
ylabel({['PC',num2str(PC),' score']})
set(gca,'FontSize',12)
t.FontSize = 12;
tString = t.String;
t.String = {[tString{1},' ',tString{2}]};
t.Position(1) = sum(xlim)/2;
t.HorizontalAlignment='center';
ylimits = ylim;
t.Position(2) = find_point_on_line(ylimits(1),ylimits(2),.95);
colormap(turbo(256))
end
end




Term=mean(PreVsTerm.TERM_Thal_TC);
Preterm=mean(PreVsTerm.PRETERM_Thal_TC);
TermN=mean(PreVsTerm.TERM_Thal_TC_norm);
PretermN=mean(PreVsTerm.PRETERM_Thal_TC_norm);
scatter(Preterm-Term,PretermN-TermN)



Term=mean(PreVsTerm.TERM_Cort_TC{3});
Preterm=mean(PreVsTerm.PRETERM_Cort_TC{3});
TermN=mean(zscore(PreVsTerm.TERM_Cort_TC_norm{3}));
PretermN=mean(zscore(PreVsTerm.PRETERM_Cort_TC_norm{3}));
scatter(Preterm-Term,PretermN-TermN)

[h,p,ci,stats] = ttest2((PreVsTerm.PRETERM_Cort_TC_norm{3}),(PreVsTerm.TERM_Cort_TC_norm{3}));


[h,p,ci,stats] = ttest2((PreVsTerm.PRETERM_Cort_TC{3}),(PreVsTerm.TERM_Cort_TC{3}));


g = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = double(g.vertices);
surface.faces = double(g.faces);

SurfData = zeros(length(medwallmask),1);

TERMvsPRE_Cort_norm_diff = mean(PreVsTerm.PRETERM_Cort_TC_norm{3})-mean(PreVsTerm.TERM_Cort_TC_norm{3});

TERMvsPRE_Cort_norm_diff = stats.tstat.*h;

SurfData(medwallmask) = TERMvsPRE_Cort_norm_diff;

plotRange = max(abs(SurfData));

SurfData(SurfData==0)=NaN;

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),'Mean difference',[-plotRange plotRange]);
