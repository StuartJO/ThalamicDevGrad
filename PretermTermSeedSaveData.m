
if Weighted == 0
    suf = '_wei';
elseif Weighted == 1
    suf = '';
end


nCortVerts = sum(medwallmask);

mkdir ./outputs/Thr100/Weighted/PvTSeed/

SeedData.Term = zeros(length(MatchedTERM),nCortVerts);
SeedData.Preterm = zeros(length(PRETERM),nCortVerts);
for s = 1:Nseed
save(['./outputs/Thr100/Weighted/PvTSeed/Seed',num2str(s),'.mat'],'-struct','SeedData');
end
%PRETERM_SIGTC = zeros(length(SigSeeds),28766,length(PRETERM));

PRETERMsum = zeros(Nseed,nCortVerts);

f = waitbar(0,'Finished 0/102');

t = 0;

nPreterm = length(PRETERM);
nTerm = length(MatchedTERM);

for i = 1:length(PRETERM)
    tic
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

for s = 1:Nseed
    SeedData = load(['./outputs/Thr100/Weighted/PvTSeed/Seed',num2str(s),'.mat']);
    SeedData.Preterm(i,:) = norm(s,:);
    save(['./outputs/Thr100/Weighted/PvTSeed/Seed',num2str(s),'.mat'],'-struct','SeedData');
end    

% PRETERM_Cort_TC(i,:) = sum(tc);
% PRETERM_Cort_TC_norm(i,:) = sum(norm);
% 
% PRETERM_Thal_TC(i,:) = sum(tc,2);
% PRETERM_Thal_TC_norm(i,:) = sum(norm,2);
    t = (toc+t);
    time_remaining = (t/i)*((nPreterm+nTerm)-i);
waitbar(i/(nPreterm+nTerm),f,['Finished ',num2str(i),'/',num2str(nPreterm+nTerm),': ',num2str(time_remaining,4),' seconds remaining'])

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

for s = 1:Nseed
    SeedData = load(['./outputs/Thr100/Weighted/PvTSeed/Seed',num2str(s),'.mat']);
    SeedData.Term(i,:) = norm(s,:);
    save(['./outputs/Thr100/Weighted/PvTSeed/Seed',num2str(s),'.mat'],'-struct','SeedData');
end    
% TERM_Cort_TC(i,:) = sum(tc);
% TERM_Cort_TC_norm(i,:) = sum(norm);
% 
% 
% TERM_Thal_TC(i,:) = sum(tc,2);
% 
% TERM_Thal_TC_norm(i,:) = sum(norm,2);

    t = (toc+t);
    time_remaining = (t./(i+Preterm))*((nPreterm+nTerm)-i-nPreterm);
waitbar((i+Preterm)/(nPreterm+nTerm),f,['Finished ',num2str(i+Preterm),'/',num2str(nPreterm+nTerm),': ',num2str(time_remaining,4),' seconds remaining'])

end
