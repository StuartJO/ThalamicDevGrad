

RescannedATterm = (thal_sub_meta.birth_age<37 & thal_sub_meta.scan_age>=37);

PRETERM = find(RescannedATterm);
TryToMatch = find(thal_sub_meta.birth_age>=37);
TryToMatch(ismember(TryToMatch,ToUse4Avg)) = [];
[~,~,SexInd] = unique(thal_sub_meta.sex);

Preterm_scan_age = thal_sub_meta.scan_age(PRETERM);
Preterm_sex = SexInd(PRETERM);
Term_scan_age = thal_sub_meta.scan_age(TryToMatch);
Term_sex = SexInd(TryToMatch);

NotMatched = true(size(Term_scan_age));
MatchedSUB = zeros(size(PRETERM));
for i = 1:length(PRETERM)
   Sex = Preterm_sex(i);
   scanAgeDiff = abs(Term_scan_age-Preterm_scan_age(i));
   scanAgeDiff(Term_sex~=Sex & ~NotMatched) = NaN;
   [~,Ind] = nanmin(scanAgeDiff);   
   MatchedSUB(i) = TryToMatch(Ind);
   NotMatched(Ind) = 0;
end

PretermIND = ones(length(PRETERM),1);
TermIND = ones(length(TryToMatch),1)*2;
GrpInd = [PretermIND;TermIND];
AgeAtScan = [Preterm_scan_age ;Term_scan_age];

% 
% for i = 1:800
%     ancova_p2(:,i) = anovan([PC1_thal_all(i,TryToMatch)';PC1_thal_all(i,Rescanned)'],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
% end
% 
% MeanDiff = mean(PC1_thal_all(:,TryToMatch)')-mean(PC1_thal_all(:,Rescanned)');
% 
% MeanDiff = mean(PC1_thal_all(:,Rescanned)')-mean(PC1_thal_all(:,TryToMatch)');
% 
% fdr_ancova = BF_FDR(ancova_p2(1,:),.05)';
% up = max(abs(MeanDiff.*fdr_ancova));
% PlotThalGradientSlices(MeanDiff.*fdr_ancova,vox_coords,cmocean('Balance'),['Term vs preterm score difference'],4,[-up up]);
% 


PretermIND = ones(length(PRETERM),1);
TermIND = ones(length(MatchedSUB),1)*2;
GrpInd = [PretermIND;TermIND];
AgeAtScan = [Preterm_scan_age ;thal_sub_meta.scan_age(MatchedSUB)];
clear ancova_p2
for i = 1:800
    ancova_p2(:,i) = anovan([PC1_thal_all(i,MatchedSUB)';PC1_thal_all(i,PRETERM)'],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
end

MeanDiff = mean(PC1_thal_all(:,PRETERM)')-mean(PC1_thal_all(:,MatchedSUB)');

%MeanDiff = mean(PC1_thal_all(:,MatchedSUB)')-mean(PC1_thal_all(:,Rescanned)');

fdr_ancova = BF_FDR(ancova_p2(1,:),.05)';
up = max(abs(MeanDiff.*fdr_ancova));
PlotThalGradientSlices(MeanDiff.*fdr_ancova,vox_coords,cmocean('Balance'),['Preterm relative to term PC1 score difference'],4,[-up up]);

%print(['TermPretermDifferences.png'],'-dpng','-r300')

NegThalDiff = (MeanDiff.*fdr_ancova)<0;
PosThalDiff = (MeanDiff.*fdr_ancova)>0;

SurfData = zeros(length(medwallmask),1);
NegThalDiffMeanConn = nanmean(NormAvg(NegThalDiff,:));
NegThalDiffMeanConn(isnan(NegThalDiffMeanConn)) = 0;
SurfData(medwallmask) = NegThalDiffMeanConn;

load('RedBluecmaps.mat')
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,neg_bar_cmap,{'Connectivity of areas with PC1 score','decreases in preterm neonates'}); 

SurfData = zeros(length(medwallmask),1);
PosThalDiffMeanConn = nanmean(NormAvg(PosThalDiff,:));
PosThalDiffMeanConn(isnan(PosThalDiffMeanConn)) = 0;
SurfData(medwallmask) = PosThalDiffMeanConn;

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,pos_bar_cmap,{'Connectivity of areas with PC1 score','increases in preterm neonates'}); 

figure
scatter(mean(PC1_thal_all(:,MatchedSUB)'),mean(PC1_thal_all(:,PRETERM)'),36,MeanDiff,'filled')
colormap(cmocean('Balance'))
up = max(abs(MeanDiff.*fdr_ancova));
caxis([-up up])
xlabel('Term PC1 scores')
ylabel('Preterm PC1 scores')
c = colorbar;
c.Label.String = {'Preterm relative to term PC1 score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)


figure
scatter(mean(PC1_thal_all(:,MatchedSUB)'),MeanDiff,36,MeanDiff,'filled')
colormap(cmocean('Balance'))
up = max(abs(MeanDiff.*fdr_ancova));
caxis([-up up])
xlabel('Term PC1 scores')
ylabel('Preterm difference')
c = colorbar;
c.Label.String = {'Preterm relative to term PC1 score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)

axis equal


figure
scatter(PC1_thal_age_corr,MeanDiff,36,MeanDiff,'filled')
colormap(cmocean('Balance'))
up = max(abs(MeanDiff.*fdr_ancova));
caxis([-up up])
xlabel('PC1 score correlation with scan age')
ylabel('Preterm difference')
c = colorbar;
c.Label.String = {'Preterm relative to term PC1 score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)
% 
% AllThalDiff = (MeanDiff.*fdr_ancova)~=0;
% 
% AllThalDiffMeanConn = nanmean(NormAvg(AllThalDiff,:));
% AllThalDiffMeanConn(isnan(AllThalDiffMeanConn)) = 0;
% 
% ExampleSurfacePlotFunction(surface,double(medwallmask),AllThalDiffMeanConn,parula(256),['Connectivity']); 
% 
