sub_ses = readtable('ALL_sub_ses.txt');

saveoutput = 1;
makefig = 1;

SUB = sub_ses.Var1;
SES = sub_ses.Var2;

Thr = 100;
load('SeedsThr.mat','ThalConnMean')
SeedThr = ThalConnMean>=100;

Nseed = sum(SeedThr);

nSub = length(SUB);

dHCP_meta = readtable('dHCP_metadata.xlsx');

for i = 1:nSub
    I1 = find(strcmp(SUB{i},dHCP_meta.participant_id));
    I2 = find(SES(i)==dHCP_meta.session_id);
    I(i) = intersect(I1,I2);   
end

thal_sub_meta = dHCP_meta(I,:);

%ToUse = (thal_sub_meta.scan_age-thal_sub_meta.birth_age)<=40;
ToUse = (thal_sub_meta.scan_age-thal_sub_meta.birth_age)<=40;

[uniqueCells, ~, duplicatesIdx] = unique(SUB, 'stable');
duplicateSubs = uniqueCells(histc(duplicatesIdx, 1:max(duplicatesIdx)) > 1);

for i = 1:length(duplicateSubs)
DupIND = find(ismember(SUB,duplicateSubs{i}));
[~,DupINDmin] = max(thal_sub_meta.scan_age(DupIND)-thal_sub_meta.birth_age(DupIND));

%[~,DupINDmin] = max(thal_sub_meta.scan_age(DupIND)-thal_sub_meta.birth_age(DupIND));

ToUse(DupIND(DupINDmin)) = false;

end

IsTerm = thal_sub_meta.birth_age>=37;

IsBestScan = thal_sub_meta.radiology_score==1;  

PossibleTrain = logical(ToUse.*IsTerm.*IsBestScan);

[~,OldestScanAgeOrd] = sort(thal_sub_meta.scan_age,'descend');

PossibleTrainOrd = PossibleTrain(OldestScanAgeOrd);

OldestScanAgeOrd(~PossibleTrainOrd) = [];

Train = OldestScanAgeOrd(1:20);

Test = find(ToUse);
Test(ismember(Test,Train)) = [];

% [OldestAge,Oldest] = sort(thal_sub_meta.birth_age,'descend');
% Oldest(ismember(Oldest,find(~ToUse))) = [];
% Train = Oldest(1:20);
% Test = Oldest(21:end);

Test_SUB = SUB(Test);
Test_SES = SES(Test);

scan_age = thal_sub_meta.scan_age(Test);
birth_age = thal_sub_meta.birth_age(Test);

%save('AgeData.mat','scan_age','birth_age')

Ntest = length(Test);

AllScore = cell(length(SUB),1);
AllCoeff = cell(length(SUB),1);

AllExpl = zeros(Ntest,Nseed-1);

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');

medwallmask = logical(medmaskgii.cdata);
Nverts = length(medwallmask);
Nverts_nomed = sum(medwallmask);

PCAOUTPUT_DIR = ['C:\Users\Stuart\Documents\GitHub\ThalamicDevGrad\outputs\Thr',num2str(Thr)];

mkdir(PCAOUTPUT_DIR)

mkdir([PCAOUTPUT_DIR,'\Weighted'])
mkdir([PCAOUTPUT_DIR,'\Unweighted'])

addpath(PCAOUTPUT_DIR)



TrainData = zeros(Nseed,length(Nverts_nomed));

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

if ismember(i,Train)
    TrainData = TrainData+data_nonmed;
    %TrainDataNorm = TrainDataNorm+Norm;
% AllExpl(i,:) = explained;    
% AllScore{i} = score(:,1:5);
% AllCoeff{i} = coeff(:,1:5);
end
disp(num2str(i))
end

TrainDataAvg = TrainData/length(Train);
save([PCAOUTPUT_DIR,'\Weighted\Avg.mat'],'TrainDataAvg')


TrainData = zeros(Nseed,length(Nverts_nomed));

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

if ismember(i,Train)
    TrainData = TrainData+data_nonmed;
    %TrainDataNorm = TrainDataNorm+Norm;
    % AllExpl(i,:) = explained;    
    % AllScore{i} = score(:,1:5);
    % AllCoeff{i} = coeff(:,1:5);
end
disp(num2str(i))
end

TrainDataAvg = TrainData/length(Train);
save([PCAOUTPUT_DIR,'\Unweighted\Avg.mat'],'TrainDataAvg')

% for i = 1:length(SUB)
% 
%     score_ = zeros(Nseed,5);
%     coeff_ = zeros(Nverts_nomed,5);
% 
%     for j = 1:5
%        score_(:,j) = dlmread([PCAOUTPUT_DIR,'\PC',num2str(j),'scores\',SUB{i},'_',num2str(SES(i)),'.txt']);
%        coeff_(:,j) = dlmread([PCAOUTPUT_DIR,'\PC',num2str(j),'coeffs\',SUB{i},'_',num2str(SES(i)),'.txt']);
%     end
%     AllExpl(i,:) = dlmread([PCAOUTPUT_DIR,'\PCvariance\',SUB{i},'_',num2str(SES(i)),'.txt']);
% 
%     AllScore{i} = score_;
%     AllCoeff{i} = coeff_;
% end



% data = dlmread(['D:/TC_connectivity/CC00063AN06_15102_thal_conn_verts_wei.txt']);
% 

%data = readmatrix(['D:/TC_connectivity/CC00063AN06_15102_thal_conn_verts_wei.txt'],'Range','1:1');
