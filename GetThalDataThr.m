Thr = 100;
load('SeedsThr.mat','ThalConnMean')
SeedThr = ThalConnMean>=Thr;

Nseed = sum(SeedThr);

AllScore = cell(length(SUB),1);
AllCoeff = cell(length(SUB),1);

AllExpl = zeros(Ntest,Nseed-1);

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');

medwallmask = logical(medmaskgii.cdata);
Nverts = length(medwallmask);
Nverts_nomed = sum(medwallmask);

PCAOUTPUT_DIR = ['.\outputs\Thr',num2str(Thr)];

mkdir(PCAOUTPUT_DIR)

mkdir([PCAOUTPUT_DIR,'\Weighted'])
mkdir([PCAOUTPUT_DIR,'\Unweighted'])

addpath(PCAOUTPUT_DIR)

TemplateData = zeros(Nseed,length(Nverts_nomed));

for i = 1:length(SUB)

data = dlmread(['D:/TC_connectivity/',SUB{i},'_',num2str(SES(i)),'_thal_conn_verts_wei.txt']);
data_nonmed = data(SeedThr,medwallmask);
data_nonmed(isnan(data_nonmed)) = 0;
Norm = BF_NormalizeMatrix(data_nonmed,'scaledSigmoid');
Norm(isnan(Norm)) = 0;
[coeff,score,~,~,varExpln] = pca(Norm);
ThalConn = sum(data_nonmed,2);
CortConn = sum(data_nonmed,1);
ThalConnNorm = sum(Norm,2);
CortConnNorm = sum(Norm,1);

coeff1_5 = coeff(:,1:5);
score1_5 = score(:,1:5);

save([PCAOUTPUT_DIR,'\Weighted\',SUB{i},'_',num2str(SES(i)),'.mat'],'coeff1_5','score1_5','CortConnNorm','ThalConnNorm','ThalConn','CortConn','varExpln')

if ismember(i,TemplateInd)
    TemplateData = TemplateData+data_nonmed;
end
disp(num2str(i))
end

TemplateDataAvg = TemplateData/length(TemplateInd);
save([PCAOUTPUT_DIR,'\Weighted\Avg.mat'],'TemplateDataAvg')


TemplateData = zeros(Nseed,length(Nverts_nomed));

for i = 1:length(SUB)

data = dlmread(['D:/TC_connectivity/',SUB{i},'_',num2str(SES(i)),'_thal_conn_verts.txt']);
data_nonmed = data(SeedThr,medwallmask);
data_nonmed(isnan(data_nonmed)) = 0;
Norm = BF_NormalizeMatrix(data_nonmed,'scaledSigmoid');
Norm(isnan(Norm)) = 0;
[coeff,score,~,~,varExpln] = pca(Norm);
ThalConn = sum(data_nonmed,2);
CortConn = sum(data_nonmed,1);
ThalConnNorm = sum(Norm,2);
CortConnNorm = sum(Norm,1);

coeff1_5 = coeff(:,1:5);
score1_5 = score(:,1:5);

save([PCAOUTPUT_DIR,'\Unweighted\',SUB{i},'_',num2str(SES(i)),'.mat'],'coeff1_5','score1_5','CortConnNorm','ThalConnNorm','ThalConn','CortConn','varExpln')

if ismember(i,TemplateInd)
    TemplateData = TemplateData+data_nonmed;
end
disp(num2str(i))
end

TemplateDataAvg = TemplateData/length(TemplateInd);
save([PCAOUTPUT_DIR,'\Unweighted\Avg.mat'],'TemplateIndDataAvg')