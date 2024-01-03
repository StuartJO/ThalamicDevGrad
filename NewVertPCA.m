sub_ses = readtable('ALL_sub_ses.txt');

saveoutput = 1;
makefig = 1;

SUB = sub_ses.Var1;
SES = sub_ses.Var2;

dHCP_meta = readtable('dHCP_metadata.xlsx');

for i = 1:363
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

save('AgeData.mat','scan_age','birth_age')

Ntest = length(Test);

AllScore = cell(length(SUB),1);
AllCoeff = cell(length(SUB),1);

AllExpl = zeros(Ntest,800);

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');

medwallmask = logical(medmaskgii.cdata);
Nverts = length(medwallmask);
Nverts_nomed = sum(medwallmask);

PCAOUTPUT_DIR = 'C:\Users\Stuart\Documents\GitHub\ThalamicDevGrad\outputs';
addpath(PCAOUTPUT_DIR)

for i = 1:length(SUB)
ThalMeanConn(i,:) = dlmread([PCAOUTPUT_DIR,'\CortPropConn\',SUB{i},'_',num2str(SES(i)),'.txt']);
end



for i = 1:length(SUB)
    
    score_ = zeros(800,5);
    coeff_ = zeros(Nverts_nomed,5);
    
    for j = 1:5
       score_(:,j) = dlmread([PCAOUTPUT_DIR,'\PC',num2str(j),'scores\',SUB{i},'_',num2str(SES(i)),'.txt']);
       coeff_(:,j) = dlmread([PCAOUTPUT_DIR,'\PC',num2str(j),'coeffs\',SUB{i},'_',num2str(SES(i)),'.txt']);
    end
    AllExpl(i,:) = dlmread([PCAOUTPUT_DIR,'\PCvariance\',SUB{i},'_',num2str(SES(i)),'.txt']);
    
    AllScore{i} = score_;
    AllCoeff{i} = coeff_;
end

TestExpl = AllExpl(Test,:);
TestScore = AllScore(Test);
TestCoeff = AllCoeff(Test);

TrainScore = dlmread('Train_avg_scores.txt');
TrainCoeff = dlmread('Train_avg_coeffs.txt');

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

TrainCoeffz_cort = zeros(Nverts,800);
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
print(['NEWVERTEX_PC1_loading_corr_scan_age_region.png'],'-dpng','-r300')
end


ExampleSurfacePlotFunction(surface,double(medwallmask),TrainCoeffz_cort(:,1),turbo(256),['Term template PC1 loading']);    

if makefig
print(['NEWVERTEX_Group_PC1_loading.png'],'-dpng','-r300')
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
print(['NEWVERTEX_Group_PC1_loading_vs_scan_age_corr.png'],'-dpng','-r300')
end

[PC1_thal_age_corr,PC1_thal_age_corr_p] = corr(thal_sub_meta.scan_age(Test),PC1_thal','Type',CorrType);
fdr_thal = BF_FDR(PC1_thal_age_corr_p,.05)';
CorrRange = max(abs(PC1_thal_age_corr(:)));
up = ceil(CorrRange * 10) / 10;

vox_coords = dlmread('thal_seed_1.75mm_vox_coords.txt');

PlotThalGradientSlices((TrainScorez(:,1)),vox_coords,turbo(256),['Term template PC1 score'],4);
if makefig
print(['NEWVERTEX_Group_PC1_score.png'],'-dpng','-r300')
end
PlotThalGradientSlices(PC1_thal_age_corr.*fdr_thal,vox_coords,cmocean('Balance'),'PC1 score correlation with scan age',4,[-up up]);
if makefig
print(['NEWVERTEX_PC1_score_corr_scan_age_region.png'],'-dpng','-r300')
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
print(['NEWVERTEX_Group_PC1_score_vs_scan_age_corr.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(thal_sub_meta.scan_age(Test),PC1_corr_wAvg);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({'Correlation with','term template (loading)'})
set(gca,'FontSize',18)
if makefig
print(['NEWVERTEX_PC1_loading_corr_scan_age.png'],'-dpng','-r300')
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
print(['NEWVERTEX_PC1_score_corr_scan_ageALPHA.png'],'-dpng','-r300')
end

figure %figure('Position',[189   318   946   688])
s = scatterfit(thal_sub_meta.scan_age(Test),TestExpl(:,1));
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({'PC1 variance explained'})
set(gca,'FontSize',18)

if saveoutput
save('NewPCAResults.mat','PC1_thal','PC1_cort','PC1_thal_age_corr','PC1_thal_age_corr_p','PC1_age_corr','PC1_age_corr_p','medwallmask','TrainScorez','TrainCoeffz')
end


% 
% figure
% plot(TrainExpl(1:6)*100,'k')
% hold on
% scatter(1:6,TrainExpl(1:6)*100,36,'filled','r')
% xlim([.5 6.5])
% xlabel('PC #')
% ylabel('% variance explained')
% set(gca,'FontSize',18)
% print(['VARIANCE_EXPLAINED.png'],'-dpng','-r300')

% FlucAmp = NaN(length(SUB),Nverts);
% for i = 1:length(SUB)
%     filename = ['sub-',SUB{i},'_ses-',num2str(SES(i)),'_hemi-left_desc-FuncAmp.shape.gii'];
%     if exist(filename,'file') == 2
%     g = gifti(filename);
%     FlucAmp(i,:) = double(g.cdata);
%     end    
% end
% 
% TrainFlucAmp = nanmean( FlucAmp(ToUse4Avg,:));
% 
% ExampleSurfacePlotFunction(surface,double(medwallmask),TrainFlucAmp,hot(256),['Term template Fluctuation Amplitude']);    
% print(['FluctuationAmplitude.png'],'-dpng','-r300')
% 
% figure
% scatterfit(TrainCoeffz(:,1),TrainFlucAmp(medwallmask))
% xlabel('PC1 loading term template')
% ylabel({'Term template','Fluctuation Amplitude'})
% set(gca,'FontSize',18)
% print(['FlucAmp_PC1_corr.png'],'-dpng','-r300')
% 
% FlucAmpTest = FlucAmp(Test,medwallmask);
% TestAge = thal_sub_meta.scan_age(Test);
% FlucAmpTestNonNaNInd = ~isnan(sum(FlucAmpTest,2));
% 
% [FlucAmp_age_corr,FlucAmp_age_corr_p] = corr(TestAge(FlucAmpTestNonNaNInd), FlucAmpTest(FlucAmpTestNonNaNInd,:),'Type',CorrType);
% 
% FlucAmp_age_corr_cort = zeros(Nverts,1);
% FlucAmp_age_corr_cort(medwallmask) = FlucAmp_age_corr;   
% ExampleSurfacePlotFunction(surface,double(medwallmask),FlucAmp_age_corr_cort,cmocean('Balance'),['Fluctuation Amplitude with age']);    
% print(['FlucAmp_age_corr.png'],'-dpng','-r300')
% 
% figure
% scatterfit(PC1_age_corr,FlucAmp_age_corr)
% xlabel('PC1 loading correlation with scan age')
% ylabel({'Fluctuation Amplitude','correlation with scan age'})
% set(gca,'FontSize',18)
% print(['FlucAmp_PC1_age_corr.png'],'-dpng','-r300')

% 
% 
% scatterfit(TrainCoeffz(:,1),PC1_cort(:,1))
% 
% s = scatterfit(vox_coords(:,2),PC1_thal_age_corr);
% s.MarkerFaceAlpha = .5;
