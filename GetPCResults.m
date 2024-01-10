PC=3;

makefig = 1;

for i = 1:length(SUB)

PC_thal_all(:,i) = zscore(aligned{i}(:,PC));
PC_aligned_all(:,i) =  zscore((coeff_aligned_all{i}(:,PC)));

end

PC_thal = PC_thal_all(:,Test);
PC_cort = PC_aligned_all(:,Test);

TrainScorez = zscore(TrainScore);
TrainCoeffz = zscore(TrainCoeff);

TrainCoeffz_cort = zeros(Nverts,799);
TrainCoeffz_cort(medwallmask,:) = TrainCoeffz;

gL = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = gL.vertices;
surface.faces = gL.faces;

CorrType = 'Pearson';

[PC_corr_wAvg,PC_corr_wAvg_p] = corr(TrainCoeffz(:,PC),PC_cort,'Type',CorrType);

[PC_corr_wAvg_thal,PC_corr_wAvg_thal_p] = corr(TrainScorez(:,PC),PC_thal,'Type',CorrType);

[PC_age_corr,PC_age_corr_p] = corr(thal_sub_meta.scan_age(Test),PC_cort','Type',CorrType);

CorrRange = max(abs(PC_age_corr(:)));
up = ceil(CorrRange * 10) / 10;

fdr = BF_FDR(PC_age_corr_p,.05)';

PC_age_corr_cort = zeros(Nverts,1);
PC_age_corr_cort(medwallmask) = PC_age_corr.*fdr;
PC_age_corr_cort(medwallmask) = PC_age_corr;   
ExampleSurfacePlotFunction(surface,double(medwallmask),PC_age_corr_cort,cmocean('Balance'),['PC loading correlation with scan age'],[-up up]);    

if makefig
print(['VERTWEI_PC',num2str(PC),'_loading_corr_scan_age_region.png'],'-dpng','-r300')
end


ExampleSurfacePlotFunction(surface,double(medwallmask),TrainCoeffz_cort(:,PC),turbo(256),['Term template PC',num2str(PC),' loading']);    

if makefig
print(['VERTWEI_Group_PC',num2str(PC),'_loading.png'],'-dpng','-r300')
end

figure %figure('Position',[189   318   946   688])
s = scatter(TrainCoeffz(:,PC),PC_age_corr,36,PC_age_corr.*fdr,'filled');
s.MarkerFaceAlpha=1;
colormap(cmocean('Balance'))
maxCorr = max(abs(PC_age_corr.*fdr));
caxis([-maxCorr maxCorr])
xlabel(['Term template PC',num2str(PC),' loading'])
ylabel({['PC',num2str(PC),' loading correlation'],'with scan age'})
set(gca,'FontSize',18)
c = colorbar;
c.Label.String = {['PC',num2str(PC),' loading correlation with scan age']};
c.Label.FontSize = 14;

if makefig
print(['VERTWEI_Group_PC',num2str(PC),'_loading_vs_scan_age_corr.png'],'-dpng','-r300')
end

[PC_thal_age_corr,PC_thal_age_corr_p] = corr(thal_sub_meta.scan_age(Test),PC_thal','Type',CorrType);
fdr_thal = BF_FDR(PC_thal_age_corr_p,.05)';
CorrRange = max(abs(PC_thal_age_corr(:)));
up = ceil(CorrRange * 10) / 10;

vox_coords = dlmread('thal_seed_1.75mm_vox_coords.txt');

PlotThalGradientSlices((TrainScorez(:,PC)),vox_coords,turbo(256),['Term template PC',num2str(PC),' score'],4);
if makefig
print(['VERTWEI_Group_PC',num2str(PC),'_score.png'],'-dpng','-r300')
end
PlotThalGradientSlices(PC_thal_age_corr.*fdr_thal,vox_coords,cmocean('Balance'),'PC score correlation with scan age',4,[-up up]);
if makefig
print(['VERTWEI_PC',num2str(PC),'_score_corr_scan_age_region.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatter(TrainScorez(:,PC),PC_thal_age_corr,36,PC_thal_age_corr.*fdr_thal,'filled');
s.MarkerFaceAlpha=.25;
colormap(cmocean('Balance'))
caxis([-up up])
xlabel(['Term template PC',num2str(PC),' score'])
ylabel({['PC',num2str(PC),' score correlation'],'with scan age'})
set(gca,'FontSize',18)
c = colorbar;
c.Label.String = {['PC',num2str(PC),' score correlation with scan age']};
c.Label.FontSize = 14;
if makefig
print(['VERTWEI_Group_PC',num2str(PC),'_score_vs_scan_age_corr.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(thal_sub_meta.scan_age(Test),PC_corr_wAvg);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({'Correlation with','term template (loading)'})
set(gca,'FontSize',18)
if makefig
print(['VERTWEI_PC',num2str(PC),'_loading_corr_scan_age.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(thal_sub_meta.scan_age(Test),PC_corr_wAvg_thal);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({'Correlation with','term template (score)'})
set(gca,'FontSize',18)

if makefig

end

s.MarkerFaceAlpha=.25;

if makefig
print(['VERTWEI_PC',num2str(PC),'_score_corr_scan_ageALPHA.png'],'-dpng','-r300')
end

% figure %figure('Position',[189   318   946   688])
% s = scatterfit(thal_sub_meta.scan_age(Test),TestExpl(:,1));
% s.MarkerFaceAlpha=1;
% xlabel('Scan age (weeks)')
% ylabel({'PC variance explained'})
% set(gca,'FontSize',18)

if saveoutput
save('VERTWEI_PCAResults.mat','PC_thal','PC_cort','PC_thal_age_corr','PC_thal_age_corr_p','PC_age_corr','PC_age_corr_p','medwallmask','TrainScorez','TrainCoeffz')
end


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
%     ancova_p2(:,i) = anovan([PC_thal_all(i,TryToMatch)';PC_thal_all(i,Rescanned)'],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
% end
% 
% MeanDiff = mean(PC_thal_all(:,TryToMatch)')-mean(PC_thal_all(:,Rescanned)');
% 
% MeanDiff = mean(PC_thal_all(:,Rescanned)')-mean(PC_thal_all(:,TryToMatch)');
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
    ancova_p2(:,i) = anovan([PC_thal_all(i,MatchedSUB)';PC_thal_all(i,PRETERM)'],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
end

MeanDiff = mean(PC_thal_all(:,PRETERM)')-mean(PC_thal_all(:,MatchedSUB)');

%MeanDiff = mean(PC_thal_all(:,MatchedSUB)')-mean(PC_thal_all(:,Rescanned)');

fdr_ancova = BF_FDR(ancova_p2(1,:),.05)';
up = max(abs(MeanDiff.*fdr_ancova));
PlotThalGradientSlices(MeanDiff.*fdr_ancova,vox_coords,cmocean('Balance'),['Preterm relative to term PC',num2str(PC),' score difference'],4,[-up up]);

%print(['TermPretermDifferences.png'],'-dpng','-r300')

NegThalDiff = (MeanDiff.*fdr_ancova)<0;
PosThalDiff = (MeanDiff.*fdr_ancova)>0;

SurfData = zeros(length(medwallmask),1);
NegThalDiffMeanConn = nanmean(Norm(NegThalDiff,:));
NegThalDiffMeanConn(isnan(NegThalDiffMeanConn)) = 0;
SurfData(medwallmask) = NegThalDiffMeanConn;

%load('RedBluecmaps.mat')

RdBlu_cmap= cmocean('Balance',1024);

neg_bar_cmap = flipud(RdBlu_cmap(1:512,:));
pos_bar_cmap = RdBlu_cmap(513:1024,:);

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,neg_bar_cmap,{['Connectivity of areas with PC',num2str(PC),' score'],'decreases in preterm neonates'}); 

SurfData = zeros(length(medwallmask),1);
PosThalDiffMeanConn = nanmean(Norm(PosThalDiff,:));
PosThalDiffMeanConn(isnan(PosThalDiffMeanConn)) = 0;
SurfData(medwallmask) = PosThalDiffMeanConn;

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,pos_bar_cmap,{['Connectivity of areas with PC',num2str(PC),' score'],'increases in preterm neonates'}); 

figure
scatter(mean(PC_thal_all(:,MatchedSUB)'),mean(PC_thal_all(:,PRETERM)'),36,MeanDiff,'filled')
colormap(cmocean('Balance'))
up = max(abs(MeanDiff.*fdr_ancova));
caxis([-up up])
xlabel(['Term PC',num2str(PC),' scores'])
ylabel(['Preterm PC',num2str(PC),' scores'])
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)


figure
scatter(mean(PC_thal_all(:,MatchedSUB)'),MeanDiff,36,MeanDiff,'filled')
colormap(cmocean('Balance'))
up = max(abs(MeanDiff.*fdr_ancova));
caxis([-up up])
xlabel(['Term PC',num2str(PC),' scores'])
ylabel('Preterm difference')
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)

axis equal

figure
scatter(PC_thal_age_corr,MeanDiff,36,MeanDiff,'filled')
colormap(cmocean('Balance'))
up = max(abs(MeanDiff.*fdr_ancova));
caxis([-up up])
xlabel(['PC',num2str(PC),' score correlation with scan age'])
ylabel('Preterm difference')
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)
