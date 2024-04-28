%% Get the meta data for the subjects to use

sub_ses = readtable('ALL_sub_ses.txt');

SUB = sub_ses.Var1;
SES = sub_ses.Var2;

ToSave.SES = SES;
ToSave.SUB = SUB;

nSub = length(SUB);

dHCP_meta = readtable('dHCP_metadata.xlsx');

% From all possible subjects, pull out the data for only those passing QC

for i = 1:nSub
    I1 = find(strcmp(SUB{i},dHCP_meta.participant_id));
    I2 = find(SES(i)==dHCP_meta.session_id);
    I(i) = intersect(I1,I2);   
end

thal_sub_meta = dHCP_meta(I,:);

%ToUse = (thal_sub_meta.scan_age-thal_sub_meta.birth_age)<=40;
ToUse = (thal_sub_meta.scan_age-thal_sub_meta.birth_age)<=40;

% Find subjects scanned more than once, use the earliest scan if so

[uniqueCells, ~, duplicatesIdx] = unique(SUB, 'stable');
duplicateSubs = uniqueCells(histc(duplicatesIdx, 1:max(duplicatesIdx)) > 1);

for i = 1:length(duplicateSubs)
    
DupIND = find(ismember(SUB,duplicateSubs{i}));
[~,DupINDmin] = max(thal_sub_meta.scan_age(DupIND)-thal_sub_meta.birth_age(DupIND));

ToUse(DupIND(DupINDmin)) = false;

end

% Get an intefer representation of sex identify. 1 = female, 2 = male
[~,~,SexInd] = unique(thal_sub_meta.sex);

ToSave.NonDup_scan_age = thal_sub_meta.scan_age(ToUse);
ToSave.NonDup_birth_age = thal_sub_meta.birth_age(ToUse);
ToSave.NonDup_sex = SexInd(ToUse);
ToSave.NonDupInd = find(ToUse);

%% Find subjects to make the term template

IsTerm = thal_sub_meta.birth_age>=37;
IsBestScan = thal_sub_meta.radiology_score==1;  
PossibleTemplate = logical(ToUse.*IsTerm.*IsBestScan);

[~,OldestScanAgeOrd] = sort(thal_sub_meta.scan_age,'descend');
PossibleTemplateOrd = PossibleTemplate(OldestScanAgeOrd);
OldestScanAgeOrd(~PossibleTemplateOrd) = [];
TemplateInd = OldestScanAgeOrd(1:20);

NonTemplateInd = find(ToUse);
NonTemplateInd(ismember(NonTemplateInd,TemplateInd)) = [];

ToSave.NonTemplate_SUB = SUB(NonTemplateInd);
ToSave.NonTemplate_SES = SES(NonTemplateInd);

ToSave.Template_SUB = SUB(TemplateInd);
ToSave.Template_SES = SES(TemplateInd);

ToSave.NonTemplate_scan_age = thal_sub_meta.scan_age(NonTemplateInd);
ToSave.NonTemplate_birth_age = thal_sub_meta.birth_age(NonTemplateInd);
ToSave.NonTemplate_sex = SexInd(NonTemplateInd);

ToSave.Template_scan_age = thal_sub_meta.scan_age(TemplateInd);
ToSave.Template_birth_age = thal_sub_meta.birth_age(TemplateInd);
ToSave.Template_sex = SexInd(TemplateInd);

ToSave.TemplateInd = TemplateInd;
ToSave.NonTemplateInd = NonTemplateInd;
%% Find preterm and matched terms for preterm analysis

% Find preterm neonates scanned at term
RescannedATterm = (thal_sub_meta.birth_age<37 & thal_sub_meta.scan_age>=37);

PretermTermScanInd = find(RescannedATterm);
% Find term neonates
TryToMatch = find(thal_sub_meta.birth_age>=37);
% Don't try to match with any in the template!
TryToMatch(ismember(TryToMatch,TemplateInd)) = [];

ToSave.PretermTermScan_SUB = SUB(PretermTermScanInd);
ToSave.PretermTermScan_SES = SES(PretermTermScanInd);

PretermTermScan_birth_age = thal_sub_meta.birth_age(PretermTermScanInd);
PretermTermScan_scan_age = thal_sub_meta.scan_age(PretermTermScanInd);
PretermTermScan_sex = SexInd(PretermTermScanInd);
Term_scan_age = thal_sub_meta.scan_age(TryToMatch);
Term_sex = SexInd(TryToMatch);

NotMatched = true(size(Term_scan_age));
MatchedTermInd = zeros(size(PretermTermScanInd));
for i = 1:length(PretermTermScanInd)
   Sex = PretermTermScan_sex(i);
   scanAgeDiff = abs(Term_scan_age-PretermTermScan_scan_age(i));
   scanAgeDiff(Term_sex~=Sex & ~NotMatched) = NaN;
   [~,Ind] = nanmin(scanAgeDiff);   
   MatchedTermInd(i) = TryToMatch(Ind);
   NotMatched(Ind) = 0;
end

ToSave.MatchedTermInd = MatchedTermInd;
ToSave.PretermTermScanInd = PretermTermScanInd;

ToSave.PretermTermScan_SUB = SUB(MatchedTermInd);
ToSave.PretermTermScan_SES = SES(MatchedTermInd);

ToSave.PretermTermScan_birth_age = PretermTermScan_birth_age;
ToSave.PretermTermScan_scan_age = PretermTermScan_scan_age;
ToSave.PretermTermScan_sex = PretermTermScan_sex;

ToSave.MatchedTerm_birth_age = thal_sub_meta.birth_age(MatchedTermInd);
ToSave.MatchedTerm_scan_age = thal_sub_meta.scan_age(MatchedTermInd);
ToSave.MatchedTerm_sex = SexInd(MatchedTermInd);

save('UsedSubData.mat','-struct','ToSave')