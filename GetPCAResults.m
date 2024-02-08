Weighted=1;

sub_ses = readtable('ALL_sub_ses.txt');

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

if Weighted == 0
WEITYPE = 'Unweighted';
WEITYPE_ABBREV='UNWEI';
elseif Weighted == 1
    WEITYPE = 'Weighted';
    WEITYPE_ABBREV='WEI';
end

Thr = 100;
load('SeedsThr.mat','ThalConnMean')
SeedThr = ThalConnMean>=100;

Nseed = sum(SeedThr);
saveoutput = 0;


Ntest = length(Test);

AllScore = cell(length(SUB),1);
AllCoeff = cell(length(SUB),1);

AllExpl = zeros(length(SUB),Nseed-1);

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');

medwallmask = logical(medmaskgii.cdata);
Nverts = length(medwallmask);
Nverts_nomed = sum(medwallmask);

%load('C:\Users\Stuart\Documents\GitHub\ThalamicDevGrad\outputs\Weighted\Avg.mat')

PCAOUTPUT_DIR = ['C:\Users\Stuart\Documents\GitHub\ThalamicDevGrad\outputs\Thr',num2str(Thr)];

load([PCAOUTPUT_DIR,'\',WEITYPE,'\Avg.mat'])

Norm = BF_NormalizeMatrix(TrainDataAvg,'scaledSigmoid');
Norm(isnan(Norm)) = 0;
[TrainCoeff,TrainScore,~,~,TrainExplained] = pca(Norm);

for i = 1:length(SUB)
    %inData = load(['C:\Users\Stuart\Documents\GitHub\ThalamicDevGrad\outputs\Weighted\',SUB{i},'_',num2str(SES(i)),'.mat']);
    inData = load([PCAOUTPUT_DIR,'\',WEITYPE,'\',SUB{i},'_',num2str(SES(i)),'.mat']);
    
    AllScore{i} = inData.score1_5;
    AllCoeff{i} = inData.coeff1_5;
    AllExpl(i,:) = inData.varExpln;
end


[aligned, xfms] = procrustes_alignment(AllScore,'reference',TrainScore(:,1:5));

for i = 1:length(SUB)
coeff_aligned_all{i} = AllCoeff{i}*xfms{i};
end

TrainScorez = zscore(TrainScore);
TrainCoeffz = zscore(TrainCoeff);

if saveoutput
save(['VERT',WEITYPE_ABBREV,'_Thr',num2str(Thr),'.mat'],'medwallmask','TrainScorez','TrainCoeffz','TrainExplained')
end

RescannedATterm = (thal_sub_meta.birth_age<37 & thal_sub_meta.scan_age>=37);

[uniqueCells, ~, duplicatesIdxP] = unique(SUB(RescannedATterm), 'stable');
duplicateSubsP = uniqueCells(histc(duplicatesIdxP, 1:max(duplicatesIdxP)) > 1);

for i = 1:length(duplicateSubsP)
DupIND = find(ismember(SUB,duplicateSubsP{i}));
[~,DupINDmin] = max(thal_sub_meta.scan_age(DupIND)-thal_sub_meta.birth_age(DupIND));
%[~,DupINDmin] = max(thal_sub_meta.scan_age(DupIND)-thal_sub_meta.birth_age(DupIND));
RescannedATterm(DupIND(DupINDmin)) = false;
end

PRETERM = find(RescannedATterm);
TryToMatch = find(thal_sub_meta.birth_age>=37);
TryToMatch(ismember(TryToMatch,Train)) = [];
[~,~,SexInd] = unique(thal_sub_meta.sex);

Preterm_scan_age = thal_sub_meta.scan_age(PRETERM);
Preterm_sex = SexInd(PRETERM);
Term_scan_age = thal_sub_meta.scan_age(TryToMatch);
Term_sex = SexInd(TryToMatch);

NotMatched = true(size(Term_scan_age));
MatchedTERM = zeros(size(PRETERM));
for i = 1:length(PRETERM)
   Sex = Preterm_sex(i);
   scanAgeDiff = abs(Term_scan_age-Preterm_scan_age(i));
   scanAgeDiff(Term_sex~=Sex & ~NotMatched) = NaN;
   [~,Ind] = nanmin(scanAgeDiff);   
   MatchedTERM(i) = TryToMatch(Ind);
   NotMatched(Ind) = 0;
end