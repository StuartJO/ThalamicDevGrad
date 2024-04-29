load('UsedSubData.mat')

if Weighted == 0
    suf = '_wei';
elseif Weighted == 1
    suf = '';
end

nCortVerts = sum(medwallmask);

PRETERM_Cort_TC = zeros(length(Preterm),nCortVerts);
PRETERM_Cort_TC_norm = zeros(length(PRETERM),nCortVerts);
TERM_Cort_TC = zeros(length(MatchedTERM),nCortVerts);
TERM_Cort_TC_norm = zeros(length(MatchedTERM),nCortVerts);

PRETERM_Thal_TC = zeros(length(PRETERM),Nseed);
PRETERM_Thal_TC_norm = zeros(length(PRETERM),Nseed);
TERM_Thal_TC = zeros(length(MatchedTERM),Nseed);
TERM_Thal_TC_norm = zeros(length(MatchedTERM),Nseed);

PRETERMsum = zeros(Nseed,nCortVerts);

for i = 1:length(PRETERM)
part_id = thal_sub_meta.participant_id(PRETERM(i));
ses_id = thal_sub_meta.session_id(PRETERM(i)); 
file = ['D:/TC_connectivity/',part_id{1},'_',num2str(ses_id),'_thal_conn_verts',suf,'.txt'];
% opts = detectImportOptions(file);
% opts.DataLines = [SigSeeds SigSeeds];
TC = readmatrix(file);
tc = TC(SeedThr,medwallmask);

tc(isnan(tc)) = 0;
PRETERMsum = PRETERMsum+tc;
norm = BF_NormalizeMatrix(tc,'scaledSigmoid');
norm(isnan(norm)) = 0;

PRETERM_Cort_TC(i,:) = sum(tc);
PRETERM_Cort_TC_norm(i,:) = sum(norm);

PRETERM_Thal_TC(i,:) = sum(tc,2);

PRETERM_Thal_TC_norm(i,:) = sum(norm,2);
end

TERMsum = zeros(Nseed,nCortVerts);
for i = 1:length(MatchedTERM)
part_id = thal_sub_meta.participant_id(MatchedTERM(i));
ses_id = thal_sub_meta.session_id(MatchedTERM(i)); 
file = ['D:/TC_connectivity/',part_id{1},'_',num2str(ses_id),'_thal_conn_verts',suf,'.txt'];

TC = readmatrix(file);
tc = TC(SeedThr,medwallmask);
tc(isnan(tc)) = 0;
TERMsum = TERMsum+tc;
norm = BF_NormalizeMatrix(tc,'scaledSigmoid');
norm(isnan(norm)) = 0;

TERM_Cort_TC(i,:) = sum(tc);
TERM_Cort_TC_norm(i,:) = sum(norm);

TERM_Thal_TC(i,:) = sum(tc,2);

TERM_Thal_TC_norm(i,:) = sum(norm,2);
end

save(['./outputs/TermVsPretermTC',WEITYPE_ABBREV,'_Thr',num2str(Thr),'.mat'],'TERM_Thal_TC_norm','PRETERM_Thal_TC_norm','TERM_Thal_TC','PRETERM_Thal_TC','TERM_Cort_TC_norm','TERM_Cort_TC','PRETERM_Cort_TC_norm','PRETERM_Cort_TC')

TERMavg = TERMsum./length(MatchedTERM);
PRETERMavg = PRETERMsum./length(MatchedTERM);
clear TERMsum PRETERMsum
save('./outputs/TermPretermAvg.mat','TERMavg','PRETERMavg')

