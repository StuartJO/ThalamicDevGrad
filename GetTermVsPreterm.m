
SigSeeds = find(fdr_ancova)';

PRETERM_SIGTC = zeros(length(SigSeeds),28766,length(PRETERM));

for i = 1:length(PRETERM)
part_id = thal_sub_meta.participant_id(PRETERM(i));
ses_id = thal_sub_meta.session_id(PRETERM(i)); 
file = ['D:/TC_connectivity/',part_id{1},'_',num2str(ses_id),'_thal_conn_verts_wei.txt'];
opts = detectImportOptions(file);
opts.DataLines = [SigSeeds SigSeeds];
tc = readmatrix(file,opts);
PRETERM_SIGTC(:,:,i) = tc(:,medwallmask);

% norm = BF_NormalizeMatrix(TrainDataAvg,'scaledSigmoid');
% norm(isnan(norm)) = 0;

end

TERM_SIGTC = zeros(length(SigSeeds),28766,length(MatchedSUB));

for i = 1:length(MatchedSUB)
part_id = thal_sub_meta.participant_id(MatchedSUB(i));
ses_id = thal_sub_meta.session_id(MatchedSUB(i)); 
file = ['D:/TC_connectivity/',part_id{1},'_',num2str(ses_id),'_thal_conn_verts_wei.txt'];
opts = detectImportOptions(file);
opts.DataLines = [SigSeeds SigSeeds];
tc = readmatrix(file,opts);
TERM_SIGTC(:,:,i) = tc(:,medwallmask);
end

TERM_SIGTC(isnan(TERM_SIGTC)) = 0;
PRETERM_SIGTC(isnan(PRETERM_SIGTC)) = 0;

TERM_SIGTC_ = reshape(permute(TERM_SIGTC,[3 1 2]),length(MatchedSUB),[]);

PRETERM_SIGTC_ = reshape(permute(PRETERM_SIGTC,[3 1 2]),length(PRETERM),[]);

[h,p,ci,stats] = ttest2(TERM_SIGTC_,PRETERM_SIGTC_);

tc_sig = reshape(h,length(SigSeeds),28766);

tc_sig(isnan(tc_sig)) = 0;

tc_t = reshape(stats.tstat,length(SigSeeds),28766);

tc_t(isnan(tc_t)) = 0;

tc_p = reshape(p,length(SigSeeds),28766);

tc_p(isnan(tc_p)) = 0;

tc_diff = mean(TERM_SIGTC,3)-mean(PRETERM_SIGTC,3);

%REV = permute(reshape(TERM_SIGTC_,length(MatchedSUB),length(SigSeeds),28766),[2 3 1]);

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = sum(tc_sig);
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,pos_bar_cmap,'Sig');

[h_,p_,ci_,stats_] = ttest2(squeeze(sum(TERM_SIGTC))',squeeze(sum(PRETERM_SIGTC))');

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = h_.*stats_.tstat;

plotRange = max(abs(SurfData));

SurfData(SurfData==0)=NaN;

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),'Sig',[-plotRange plotRange]);

