Weighted=1;
saveoutput = 1;

load('UsedSubData.mat')

if Weighted == 0
WEITYPE = 'Unweighted';
WEITYPE_ABBREV='UNWEI';
elseif Weighted == 1
    WEITYPE = 'Weighted';
    WEITYPE_ABBREV='WEI';
end

Thr = 100;
load('SeedsThr.mat','ThalConnMean')
SeedThr = ThalConnMean>=Thr;
Nseed = sum(SeedThr);

AllScore = cell(length(SUB),1);
AllCoeff = cell(length(SUB),1);

AllExpl = zeros(length(SUB),Nseed-1);

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');

medwallmask = logical(medmaskgii.cdata);
Nverts = length(medwallmask);
Nverts_nomed = sum(medwallmask);

%load('C:\Users\Stuart\Documents\GitHub\ThalamicDevGrad\outputs\Weighted\Avg.mat')

PCAOUTPUT_DIR = ['.\outputs\Thr',num2str(Thr)];

load([PCAOUTPUT_DIR,'\',WEITYPE,'\Avg.mat'])

Norm = BF_NormalizeMatrix(TemplateDataAvg,'scaledSigmoid');
Norm(isnan(Norm)) = 0;
[TemplateCoeff,TemplateScore,~,~,TemplateExplained] = pca(Norm);

for i = 1:length(SUB)
    
    inData = load([PCAOUTPUT_DIR,'\',WEITYPE,'\',SUB{i},'_',num2str(SES(i)),'.mat']);
    
    AllScore{i} = inData.score1_5;
    AllCoeff{i} = inData.coeff1_5;
    AllExpl(i,:) = inData.varExpln;
end

[aligned, xfms] = procrustes_alignment(AllScore,'reference',TemplateScore(:,1:5));

for i = 1:length(SUB)
coeff_aligned_all{i} = AllCoeff{i}*xfms{i};
end

TemplateScorez = zscore(TemplateScore);
TemplateCoeffz = zscore(TemplateCoeff);

if saveoutput
save(['.\outputs\VERT',WEITYPE_ABBREV,'_Thr',num2str(Thr),'.mat'],'medwallmask','TemplateScorez','TemplateCoeffz','TemplateExplained')
end
