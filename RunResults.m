sub_ses = readtable('ALL_sub_ses.txt');

saveoutput = 1;
makefig = 1;

SUB = sub_ses.Var1;
SES = sub_ses.Var2;

nSub = length(SUB);

dHCP_meta = readtable('dHCP_metadata.xlsx');

for i = 1:nSub
    I1 = find(strcmp(SUB{i},dHCP_meta.participant_id));
    I2 = find(SES(i)==dHCP_meta.session_id);
    I(i) = intersect(I1,I2);   
end

thal_sub_meta = dHCP_meta(I,:);

ToUse = (thal_sub_meta.scan_age-thal_sub_meta.birth_age)<=4;

[uniqueCells, ~, duplicatesIdx] = unique(SUB, 'stable');
duplicateSubs = uniqueCells(histc(duplicatesIdx, 1:max(duplicatesIdx)) > 1);

for i = 1:length(duplicateSubs)
DupIND = find(ismember(SUB,duplicateSubs{i}));
[~,DupINDmin] = max(thal_sub_meta.scan_age(DupIND)-thal_sub_meta.birth_age(DupIND));

ToUse(DupIND(DupINDmin)) = false;

end

[OldestAge,Oldest] = sort(thal_sub_meta.birth_age,'descend');
Oldest(ismember(Oldest,find(~ToUse))) = [];
ToUse4Avg = Oldest(1:20);
Test = Oldest(21:end);

Test_SUB = SUB(Test);
Test_SES = SES(Test);

scan_age = thal_sub_meta.scan_age(Test);
birth_age = thal_sub_meta.birth_age(Test);

Ntest = length(Test);

AllScore = cell(length(SUB),1);
AllCoeff = cell(length(SUB),1);

AllExpl = zeros(Ntest,799);

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');

medwallmask = logical(medmaskgii.cdata);
Nverts = length(medwallmask);
Nverts_nomed = sum(medwallmask);

load('C:\Users\Stuart\Documents\GitHub\ThalamicDevGrad\outputs\Weighted\Avg.mat')
Norm = BF_NormalizeMatrix(TrainDataAvg,'scaledsigmoid');
Norm(isnan(Norm)) = 0;
[TrainCoeff,TrainScore,~,~,TrainExplained] = pca(Norm);

for i = 1:length(SUB)
    inData = load(['C:\Users\Stuart\Documents\GitHub\ThalamicDevGrad\outputs\Weighted\',SUB{i},'_',num2str(SES(i)),'.mat']);
    AllScore{i} = inData.score1_5;
    AllCoeff{i} = inData.coeff1_5;
end

[aligned, xfms] = procrustes_alignment(AllScore,'reference',TrainScore(:,1:5));

for i = 1:length(SUB)
coeff_aligned_all{i} = AllCoeff{i}*xfms{i};

PC1_thal_all(:,i) = zscore(aligned{i}(:,1));
PC1_aligned_all(:,i) =  zscore((coeff_aligned_all{i}(:,1)));
PC2_thal_all(:,i) =  zscore(aligned{i}(:,2));
PC2_aligned_all(:,i) =  zscore((coeff_aligned_all{i}(:,2)));

end

PC1_thal = PC1_thal_all(:,Test);
PC1_cort = PC1_aligned_all(:,Test);

TrainScorez = zscore(TrainScore);
TrainCoeffz = zscore(TrainCoeff);

TrainCoeffz_cort = zeros(Nverts,799);
TrainCoeffz_cort(medwallmask,:) = TrainCoeffz;

gL = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = gL.vertices;
surface.faces = gL.faces;

CorrType = 'Pearson';

[PC1_corr_wAvg,PC1_corr_wAvg_p] = corr(TrainCoeffz(:,1),PC1_cort,'Type',CorrType);

[PC1_corr_wAvg_thal,PC1_corr_wAvg_thal_p] = corr(TrainScorez(:,1),PC1_thal,'Type',CorrType);

[PC1_age_corr,PC1_age_corr_p] = corr(thal_sub_meta.scan_age(Test),PC1_cort','Type',CorrType);

CorrRange = max(abs(PC1_age_corr(:)));
up = ceil(CorrRange * 10) / 10;

fdr = BF_FDR(PC1_age_corr_p,.05)';

PC1_age_corr_cort = zeros(Nverts,1);
PC1_age_corr_cort(medwallmask) = PC1_age_corr.*fdr;
PC1_age_corr_cort(medwallmask) = PC1_age_corr;   
ExampleSurfacePlotFunction(surface,double(medwallmask),PC1_age_corr_cort,cmocean('Balance'),['PC1 loading correlation with scan age'],[-up up]);    

if makefig
print(['VERTWEI_PC1_loading_corr_scan_age_region.png'],'-dpng','-r300')
end


ExampleSurfacePlotFunction(surface,double(medwallmask),TrainCoeffz_cort(:,1),turbo(256),['Term template PC1 loading']);    

if makefig
print(['VERTWEI_Group_PC1_loading.png'],'-dpng','-r300')
end

figure %figure('Position',[189   318   946   688])
s = scatter(TrainCoeffz(:,1),PC1_age_corr,36,PC1_age_corr.*fdr,'filled');
s.MarkerFaceAlpha=1;
colormap(cmocean('Balance'))
maxCorr = max(abs(PC1_age_corr.*fdr));
caxis([-maxCorr maxCorr])
xlabel('Term template PC1 loading')
ylabel({'PC1 loading correlation','with scan age'})
set(gca,'FontSize',18)
c = colorbar;
c.Label.String = {'PC1 loading correlation with scan age'};
c.Label.FontSize = 14;

if makefig
print(['VERTWEI_Group_PC1_loading_vs_scan_age_corr.png'],'-dpng','-r300')
end

[PC1_thal_age_corr,PC1_thal_age_corr_p] = corr(thal_sub_meta.scan_age(Test),PC1_thal','Type',CorrType);
fdr_thal = BF_FDR(PC1_thal_age_corr_p,.05)';
CorrRange = max(abs(PC1_thal_age_corr(:)));
up = ceil(CorrRange * 10) / 10;

vox_coords = dlmread('thal_seed_1.75mm_vox_coords.txt');

PlotThalGradientSlices((TrainScorez(:,1)),vox_coords,turbo(256),['Term template PC1 score'],4);
if makefig
print(['VERTWEI_Group_PC1_score.png'],'-dpng','-r300')
end
PlotThalGradientSlices(PC1_thal_age_corr.*fdr_thal,vox_coords,cmocean('Balance'),'PC1 score correlation with scan age',4,[-up up]);
if makefig
print(['VERTWEI_PC1_score_corr_scan_age_region.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatter(TrainScorez(:,1),PC1_thal_age_corr,36,PC1_thal_age_corr.*fdr_thal,'filled');
s.MarkerFaceAlpha=.25;
colormap(cmocean('Balance'))
caxis([-up up])
xlabel('Term template PC1 score')
ylabel({'PC1 score correlation','with scan age'})
set(gca,'FontSize',18)
c = colorbar;
c.Label.String = {'PC1 score correlation with scan age'};
c.Label.FontSize = 14;
if makefig
print(['VERTWEI_Group_PC1_score_vs_scan_age_corr.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(thal_sub_meta.scan_age(Test),PC1_corr_wAvg);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({'Correlation with','term template (loading)'})
set(gca,'FontSize',18)
if makefig
print(['VERTWEI_PC1_loading_corr_scan_age.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(thal_sub_meta.scan_age(Test),PC1_corr_wAvg_thal);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({'Correlation with','term template (score)'})
set(gca,'FontSize',18)

if makefig

end

s.MarkerFaceAlpha=.25;

if makefig
print(['VERTWEI_PC1_score_corr_scan_ageALPHA.png'],'-dpng','-r300')
end

% figure %figure('Position',[189   318   946   688])
% s = scatterfit(thal_sub_meta.scan_age(Test),TestExpl(:,1));
% s.MarkerFaceAlpha=1;
% xlabel('Scan age (weeks)')
% ylabel({'PC1 variance explained'})
% set(gca,'FontSize',18)

if saveoutput
save('VERTWEI_PCAResults.mat','PC1_thal','PC1_cort','PC1_thal_age_corr','PC1_thal_age_corr_p','PC1_age_corr','PC1_age_corr_p','medwallmask','TrainScorez','TrainCoeffz')
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
NegThalDiffMeanConn = nanmean(Norm(NegThalDiff,:));
NegThalDiffMeanConn(isnan(NegThalDiffMeanConn)) = 0;
SurfData(medwallmask) = NegThalDiffMeanConn;

%load('RedBluecmaps.mat')

RdBlu_cmap= cmocean('Balance',1024);

neg_bar_cmap = flipud(RdBlu_cmap(1:512,:));
pos_bar_cmap = RdBlu_cmap(513:1024,:);

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,neg_bar_cmap,{'Connectivity of areas with PC1 score','decreases in preterm neonates'}); 

SurfData = zeros(length(medwallmask),1);
PosThalDiffMeanConn = nanmean(Norm(PosThalDiff,:));
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
