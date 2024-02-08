

makefig = 0;
saveoutput=0;

Weighted=1;

if Weighted == 0
WEITYPE = 'Unweighted';
WEITYPE_ABBREV='_UNWEI';
elseif Weighted == 1
    WEITYPE = 'Weighted';
    WEITYPE_ABBREV='_WEI';
end

PCAOUTPUTFIG_DIR = ['C:/Users/Stuart/Documents/GitHub/ThalamicDevGrad/figures/Thr',num2str(Thr),'/',WEITYPE];

mkdir(PCAOUTPUTFIG_DIR)

for PC=1

for i = 1:length(SUB)

PC_thal_all(:,i) = zscore(aligned{i}(:,PC));
PC_cort_all(:,i) =  zscore((coeff_aligned_all{i}(:,PC)));

end

PC_thal = PC_thal_all(:,Test);
PC_cort = PC_cort_all(:,Test);

TrainCoeffz_cort = zeros(Nverts,Nseed-1);
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
ExampleSurfacePlotFunction(surface,double(medwallmask),PC_age_corr_cort,cmocean('Balance'),['PC',num2str(PC),' loading correlation with scan age'],[-up up]);    

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_loading_corr_scan_age_region.png'],'-dpng','-r300')
end


ExampleSurfacePlotFunction(surface,double(medwallmask),TrainCoeffz_cort(:,PC),turbo(256),['Term template PC',num2str(PC),' loading']);    

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_Group_PC',num2str(PC),'_loading.png'],'-dpng','-r300')
end

figure %figure('Position',[189   318   946   688])
s = scatterfit(TrainCoeffz(:,PC),PC_age_corr,36,PC_age_corr.*fdr,[],2);
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
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_Group_PC',num2str(PC),'_loading_vs_scan_age_corr.png'],'-dpng','-r300')
end

[PC_thal_age_corr,PC_thal_age_corr_p] = corr(thal_sub_meta.scan_age(Test),PC_thal','Type',CorrType);
fdr_thal = BF_FDR(PC_thal_age_corr_p,.05)';
CorrRange = max(abs(PC_thal_age_corr(:)));
up = ceil(CorrRange * 10) / 10;

vox_coords_full = dlmread('thal_seed_1.75mm_vox_coords.txt');

vox_coords = vox_coords_full(SeedThr,:);

PlotThalGradientSlices((TrainScorez(:,PC)),vox_coords,turbo(256),['Term template PC',num2str(PC),' score'],4);
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_Group_PC',num2str(PC),'_score.png'],'-dpng','-r300')
end
PlotThalGradientSlices(PC_thal_age_corr.*fdr_thal,vox_coords,cmocean('Balance'),['PC',num2str(PC),' score correlation with scan age'],4,[-up up]);
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_score_corr_scan_age_region.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(TrainScorez(:,PC),PC_thal_age_corr,36,PC_thal_age_corr.*fdr_thal,[],2);
s.MarkerFaceAlpha=1;
colormap(cmocean('Balance'))
caxis([-up up])
xlabel(['Term template PC',num2str(PC),' score'])
ylabel({['PC',num2str(PC),' score correlation'],'with scan age'})
set(gca,'FontSize',18)
c = colorbar;
c.Label.String = {['PC',num2str(PC),' score correlation with scan age']};
c.Label.FontSize = 14;
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_Group_PC',num2str(PC),'_score_vs_scan_age_corr.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(thal_sub_meta.scan_age(Test),PC_corr_wAvg);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({['Correlation with PC',num2str(PC)],'term template (loading)'})
set(gca,'FontSize',18)
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_loading_corr_scan_age.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(thal_sub_meta.scan_age(Test),PC_corr_wAvg_thal);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({['Correlation with PC',num2str(PC)],'term template (score)'})
set(gca,'FontSize',18)



s.MarkerFaceAlpha=1;

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_score_corr_scan_ageALPHA.png'],'-dpng','-r300')
end

figure('Position',[62   300   1300   335])
for i = 1:3
    subplot(1,3,i)
[s,t] = scatterfit(vox_coords(:,i),PC_thal_age_corr,36,TrainScorez(:,PC));
s.MarkerFaceAlpha=1;
cardDir={'Medial-Lateral','Anterior-posterior','Ventral-Dorsal'};
xlabel(cardDir{i})
ylabel({['PC',num2str(PC),' score correlation'],'with scan age'})
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
c = colorbar;
c.Position = [0.927282051527806	0.146778037162555	0.0151076923076924	0.778042971020951];
c.Label.String = ['Term template PC',num2str(PC),' loading'];

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_score_CardDirCorr.png'],'-dpng','-r300')
end

% figure %figure('Position',[189   318   946   688])
% s = scatterfit(thal_sub_meta.scan_age(Test),TestExpl(:,1));
% s.MarkerFaceAlpha=1;
% xlabel('Scan age (weeks)')
% ylabel({'PC variance explained'})
% set(gca,'FontSize',18)


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
TermIND = ones(length(MatchedTERM),1)*2;
GrpInd = [PretermIND;TermIND];
AgeAtScan = [Preterm_scan_age ;thal_sub_meta.scan_age(MatchedTERM)];
clear ancova_thal
for i = 1:Nseed
    ancova_thal(:,i) = anovan([PC_thal_all(i,PRETERM)';PC_thal_all(i,MatchedTERM)'],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
end

PvT_thalDiff = mean(PC_thal_all(:,PRETERM)')-mean(PC_thal_all(:,MatchedTERM)');

%MeanDiff = mean(PC_thal_all(:,MatchedSUB)')-mean(PC_thal_all(:,Rescanned)');

fdr_thal = BF_FDR(ancova_thal(1,:),.05)';


ancova_cort = zeros(Nverts_nomed,1);
F = zeros(Nverts_nomed,1);
for i = 1:Nverts_nomed
    [p,tbl,stats] = anovan([PC_cort_all(i,PRETERM)';PC_cort_all(i,MatchedTERM)'],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
ancova_cort(i,:) = p(1);
F(i,:) = tbl{2,6};
end

fdr_cort = BF_FDR(ancova_cort,.05)';

PvT_cortDiff = mean(PC_cort_all(:,PRETERM)')-mean(PC_cort_all(:,MatchedTERM)');

diff_sign = sign(PvT_cortDiff)';


if saveoutput
save(['VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_Thr',num2str(Thr),'_Results.mat'],'PC_thal','PC_cort','PC_thal_age_corr','PC_thal_age_corr_p','PC_age_corr','PC_age_corr_p','fdr_thal','ancova_thal','PvT_thalDiff','PvT_cortDiff','ancova_cort','fdr_cort')
end

up = max(abs(PvT_thalDiff.*fdr_thal));
PlotThalGradientSlices(PvT_thalDiff.*fdr_thal,vox_coords,cmocean('Balance'),['Preterm relative to term PC',num2str(PC),' score difference'],4,[-up up]);

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_thalDiff_corrected.png'],'-dpng','-r300')
end

uncorrected = ancova_thal(1,:)<.05;
up = max(abs(PvT_thalDiff.*uncorrected));
PlotThalGradientSlices(PvT_thalDiff,vox_coords,cmocean('Balance'),['Preterm relative to term PC',num2str(PC),' score difference'],4,[-up up]);

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_thalDiff.png'],'-dpng','-r300')
end

%%
val = PvT_cortDiff;
SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = val;
up = max(abs(val));
RdBlu_cmap= cmocean('Balance',1024);

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),['Preterm vs term PC',num2str(PC),' differences'],[-up up]); 

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_cortDiff.png'],'-dpng','-r300')
end

val(~fdr_cort) = NaN;
up = nanmax(abs(val));
SurfData(medwallmask) = val;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),['Preterm vs term PC',num2str(PC),' sig. differences'],[-up up]); 

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_cortDiff_corrected.png'],'-dpng','-r300')
end

%%
figure
scatter(mean(PC_cort_all(:,MatchedTERM)'),mean(PC_cort_all(:,PRETERM)'),10,PvT_cortDiff.*fdr_cort,'filled')
colormap(cmocean('Balance'))
up = max(abs(PvT_cortDiff.*fdr_cort));
caxis([-up up])
xlabel(['Term PC',num2str(PC),' loadings'])
ylabel(['Preterm PC',num2str(PC),' loadings'])
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' loadings difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_cortScatter.png'],'-dpng','-r300')
end
%%

%print(['TermPretermDifferences.png'],'-dpng','-r300')

NegThalDiff = (PvT_thalDiff.*fdr_thal)<0;
PosThalDiff = (PvT_thalDiff.*fdr_thal)>0;

SurfData = zeros(length(medwallmask),1);
NegThalDiffMeanConn = nanmean(Norm(NegThalDiff,:));
NegThalDiffMeanConn(isnan(NegThalDiffMeanConn)) = 0;
SurfData(medwallmask) = NegThalDiffMeanConn;

%load('RedBluecmaps.mat')

RdBlu_cmap= cmocean('Balance',1024);

neg_bar_cmap = flipud(RdBlu_cmap(1:512,:));
pos_bar_cmap = RdBlu_cmap(513:1024,:);

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,neg_bar_cmap,{['Connectivity of areas with PC',num2str(PC),' score'],'decreases in preterm neonates'}); 

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_cortNeg.png'],'-dpng','-r300')
end

SurfData = zeros(length(medwallmask),1);
PosThalDiffMeanConn = nanmean(Norm(PosThalDiff,:));
PosThalDiffMeanConn(isnan(PosThalDiffMeanConn)) = 0;
SurfData(medwallmask) = PosThalDiffMeanConn;

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,pos_bar_cmap,{['Connectivity of areas with PC',num2str(PC),' score'],'increases in preterm neonates'}); 

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_cortPos.png'],'-dpng','-r300')
end


figure
scatter(mean(PC_thal_all(:,MatchedTERM)'),mean(PC_thal_all(:,PRETERM)'),10,PvT_thalDiff.*fdr_ancova,'filled')
colormap(cmocean('Balance'))
up = max(abs(MeanDiff.*fdr_ancova));
caxis([-up up])
xlabel(['Term PC',num2str(PC),' scores'])
ylabel(['Preterm PC',num2str(PC),' scores'])
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_thalScatter.png'],'-dpng','-r300')
end

figure
scatter(mean(PC_thal_all(:,MatchedTERM)'),PvT_thalDiff,36,PvT_thalDiff,'filled')
colormap(cmocean('Balance'))
% up = max(abs(MeanDiff.*fdr_ancova));
% caxis([-up up])
xlabel(['Term PC',num2str(PC),' scores'])
ylabel('Preterm difference')
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)

axis equal

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSter_thalDiffm.png'],'-dpng','-r300')
end

figure
scatter(PC_thal_age_corr,PvT_thalDiff,36,PvT_thalDiff,'filled')
colormap(cmocean('Balance'))
% up = max(abs(MeanDiff.*fdr_ancova));
% caxis([-up up])
xlabel(['PC',num2str(PC),' score correlation with scan age'])
ylabel('Preterm difference')
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_thalAge.png'],'-dpng','-r300')
end

figure('Position',[62   300   1300   335])
for i = 1:3
    subplot(1,3,i)
[s,t] = scatterfit(vox_coords(:,i),TrainScorez(:,PC),36,lines(1),[],2);
s.MarkerFaceAlpha=.5;
cardDir={'Medial-Lateral','Anterior-posterior','Ventral-Dorsal'};
xlabel(cardDir{i})
ylabel(['PC',num2str(PC),' score'])
set(gca,'FontSize',12)
% t.FontSize = 12;
% tString = t.String;
% t.String = {[tString{1},' ',tString{2}]};
% t.Position(1) = sum(xlim)/2;
% t.HorizontalAlignment='center';
% ylimits = ylim;
% t.Position(2) = find_point_on_line(ylimits(1),ylimits(2),.95);
% colormap(turbo(256))
end
% c = colorbar;
% c.Position = [0.927282051527806	0.146778037162555	0.0151076923076924	0.778042971020951];
% c.Label.String = ['Term template PC',num2str(PC),' loading'];

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_score_CardDir.png'],'-dpng','-r300')
end

close all

end