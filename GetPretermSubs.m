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