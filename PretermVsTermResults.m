
if Weighted == 0
    suf = '_wei';
elseif Weighted == 1
    suf = '';
end

for PC = 1:2
PCin{PC} = load(['VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_Thr',num2str(Thr),'_Results.mat']);
end
run = 0;
if run
nCortVerts = sum(medwallmask);

PRETERM_Cort_TC_PC = cell(2,4);
PRETERM_Cort_TC_PC_norm = cell(2,4);
TERM_Cort_TC_PC = cell(2,4);
TERM_Cort_TC_PC_norm = cell(2,4);

PRETERM_Cort_TC = zeros(length(PRETERM),nCortVerts);
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

save(['TermVsPretermTC',WEITYPE_ABBREV,'_Thr',num2str(Thr),'.mat'],'TERM_Thal_TC_norm','PRETERM_Thal_TC_norm','TERM_Thal_TC','PRETERM_Thal_TC','TERM_Cort_TC_norm','TERM_Cort_TC','PRETERM_Cort_TC_norm','PRETERM_Cort_TC')

TERMavg = TERMsum./length(MatchedTERM);
PRETERMavg = PRETERMsum./length(MatchedTERM);
clear TERMsum PRETERMsum
save('TermPretermAvg.mat','TERMavg','PRETERMavg')
else
    load('TermVsPretermTC_WEI_Thr100.mat')
end

nCortVerts = sum(medwallmask);

PretermIND = ones(length(PRETERM),1);
TermIND = ones(length(MatchedTERM),1)*2;
GrpInd = [PretermIND;TermIND];
AgeAtScan = [Preterm_scan_age ;thal_sub_meta.scan_age(MatchedTERM)];

TERMvsPRE_Cort_norm = zeros(nCortVerts,3);
TERMvsPRE_Cort = zeros(nCortVerts,3);

preterm_data = PRETERM_Cort_TC_norm;
term_data =TERM_Cort_TC_norm;
for i = 1:nCortVerts
    [p,tbl] = anovan([preterm_data(:,i);term_data(:,i)],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
    TERMvsPRE_Cort_norm(i,2) = p(1);
    TERMvsPRE_Cort_norm(i,1) = tbl{2,6};
end
TERMvsPRE_Cort_norm(:,3) = mean(preterm_data)-mean(term_data);

preterm_data = PRETERM_Cort_TC;
term_data =TERM_Cort_TC;
for i = 1:nCortVerts
    [p,tbl] = anovan([preterm_data(:,i);term_data(:,i)],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
    TERMvsPRE_Cort(i,2) = p(1);
    TERMvsPRE_Cort(i,1) = tbl{2,6};
end
TERMvsPRE_Cort(:,3) = mean(preterm_data)-mean(term_data);

TERMvsPRE_Thal_norm = zeros(Nseed,3);
for i = 1:Nseed
    [p,tbl] = anovan([PRETERM_Thal_TC_norm(:,i);TERM_Thal_TC_norm(:,i)],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
    TERMvsPRE_Thal_norm(i,2) = p(1);
    TERMvsPRE_Thal_norm(i,1) = tbl{2,6};
end
TERMvsPRE_Thal_norm(:,3) = mean(PRETERM_Thal_TC_norm)-mean(TERM_Thal_TC_norm);

TERMvsPRE_Thal = zeros(Nseed,3);
for i = 1:Nseed
    [p,tbl] = anovan([PRETERM_Thal_TC(:,i);TERM_Thal_TC(:,i)],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
    TERMvsPRE_Thal(i,2) = p(1);
    TERMvsPRE_Thal(i,1) = tbl{2,6};
end
TERMvsPRE_Thal(:,3) = mean(PRETERM_Thal_TC)-mean(TERM_Thal_TC);

%%

fdr = BF_FDR(TERMvsPRE_Cort_norm(:,2),.05);
SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort_norm(:,1).*fdr.*sign(TERMvsPRE_Cort_norm(:,3));
plotRange = max(abs(SurfData));
SurfData(SurfData==0)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"Preterm-term {\itF} statistic",[-plotRange plotRange],'none');

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_cort_F_FDR.png'],'-dpng','-r300')

thal_sig = BF_FDR(TERMvsPRE_Thal_norm(:,2),.05);
%thal_sig = TERMvsPRE_Thal_norm(:,2)<.05;
ThalData2Plot = TERMvsPRE_Thal_norm(:,1).*thal_sig.*sign(TERMvsPRE_Thal_norm(:,3));
up = max(abs(ThalData2Plot));
ThalData2Plot(~thal_sig)=NaN;
PlotThalGradientSlices(ThalData2Plot,vox_coords,cmocean('Balance'),"Preterm-term {\itF} statistic",4,[-up up]);

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_thal_F_FDR.png'],'-dpng','-r300')

%%
SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort(:,1).*sign(TERMvsPRE_Cort(:,3));
plotRange = max(abs(SurfData));
SurfData(SurfData==0)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"Preterm-term {\itF} statistic",[-plotRange plotRange],'none');

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_cortRaw_F.png'],'-dpng','-r300')

ThalData2Plot = TERMvsPRE_Thal(:,1).*sign(TERMvsPRE_Thal(:,3));
up = max(abs(ThalData2Plot));
PlotThalGradientSlices(ThalData2Plot,vox_coords,cmocean('Balance'),"Preterm-term {\itF} statistic",4,[-up up]);

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_thalRaw_F.png'],'-dpng','-r300')

%%

fdr = BF_FDR(TERMvsPRE_Cort(:,2),.05);
SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort(:,1).*fdr.*sign(TERMvsPRE_Cort(:,3));
plotRange = max(abs(SurfData));
SurfData(SurfData==0)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"Preterm-term {\itF} statistic",[-plotRange plotRange],'none');

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_cortRaw_F_FDR.png'],'-dpng','-r300')

thal_sig = BF_FDR(TERMvsPRE_Thal(:,2),.05);
%thal_sig = TERMvsPRE_Thal_norm(:,2)<.05;
ThalData2Plot = TERMvsPRE_Thal(:,1).*thal_sig.*sign(TERMvsPRE_Thal(:,3));
up = max(abs(ThalData2Plot));
ThalData2Plot(~thal_sig)=NaN;
PlotThalGradientSlices(ThalData2Plot,vox_coords,cmocean('Balance'),"Preterm-term {\itF} statistic",4,[-up up]);

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_thalRaw_F_FDR.png'],'-dpng','-r300')

%%
SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort_norm(:,1).*sign(TERMvsPRE_Cort_norm(:,3));
plotRange = max(abs(SurfData));
SurfData(SurfData==0)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"Preterm-term {\itF} statistic",[-plotRange plotRange],'none');

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_cort_F.png'],'-dpng','-r300')

ThalData2Plot = TERMvsPRE_Thal_norm(:,1).*sign(TERMvsPRE_Thal_norm(:,3));
up = max(abs(ThalData2Plot));
PlotThalGradientSlices(ThalData2Plot,vox_coords,cmocean('Balance'),"Preterm-term {\itF} statistic",4,[-up up]);

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_thal_F.png'],'-dpng','-r300')

%%

figure('Position',[62   300   1300   335])
for i = 1:3
    subplot(1,3,i)
[s,t] = scatterfit(vox_coords(:,i),TERMvsPRE_Thal_norm(:,3),36,lines(1),[],2);
s.MarkerFaceAlpha=.5;
cardDir={'Medial-Lateral','Anterior-posterior','Ventral-Dorsal'};
xlabel(cardDir{i})
ylabel(['Preterm-term difference'])
set(gca,'FontSize',12)
end
% c = colorbar;
% c.Position = [0.927282051527806	0.146778037162555	0.0151076923076924	0.778042971020951];
% c.Label.String = ['Term template PC',num2str(PC),' loading'];

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'PvT_diff_CardDir.png'],'-dpng','-r300')
end



PC1 = load('VERT_WEI_PC1_Thr100_Results.mat');
PC2 = load('VERT_WEI_PC2_Thr100_Results.mat');
PC3 = load('VERT_WEI_PC3_Thr100_Results.mat');

DATA(:,1) = TrainCoeffz(:,1);
DATA(:,2) = TrainCoeffz(:,2);
DATA(:,3) = TrainCoeffz(:,3);

DATA(:,4) = PC1.PC_age_corr;
DATA(:,5) = PC2.PC_age_corr;
DATA(:,6) = PC3.PC_age_corr;

DATA_LABEL = {'PC1 template loading','PC2 template loading','PC3 template loading','PC1 age loading corr.','PC2 age loading corr.','PC3 age loading corr.'};
DATA2_LABEL = {'PC1 template score','PC2 template score','PC3 template score','PC1 age score corr.','PC2 age score corr.','PC3 age score corr.'};

DATA2(:,1) = TrainScorez(:,1);
DATA2(:,2) = TrainScorez(:,2);
DATA2(:,3) = TrainScorez(:,3);

DATA2(:,4) = PC1.PC_thal_age_corr;
DATA2(:,5) = PC2.PC_thal_age_corr;
DATA2(:,6) = PC3.PC_thal_age_corr;

SurfData = NaN(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort_norm(:,1).*sign(TERMvsPRE_Cort_norm(:,3));
figure('Position',[1 49 1382 747])
for i = 1:6
InDATA = NaN(length(medwallmask),1);
InDATA(medwallmask) = DATA(:,i);
p_perm = perm_sphere_p(InDATA,SurfData,spins','Pearson');
%figure('Position', [291 391 560 420])
subplot(2,3,i)  
if p_perm == 0
    pval_format = '{\itp_{spin}} < .001 ';   
else
    pval_format = ['{\itp_{spin}} = ', num2str(round(p_perm, 4))];
end
s = scatterfit(InDATA, SurfData, 40,lines(1),pval_format,2);
%title(['{\itr} = ', num2str(neuromap_corrs.corr(SigInd), 3), ', ', pval_format, ''], 'FontWeight', 'normal');
% xlabel(xlabel_var)
% ylabel(neuromap_corrs.name{SigInd})
% colormap(turbo)
xlabel(DATA_LABEL{i})
ylabel({"Cortical preterm-term","{\itF} statistic"})

set(gca, 'FontSize', 16)
s.MarkerEdgeColor = [0 0 0];

end

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_cort_diff_scatts.png'],'-dpng','-r300')

PREvsTERM_diff_thal = TERMvsPRE_Thal(:,1).*sign(TERMvsPRE_Thal(:,3));

% [~,~,~,stats] = ttest2(PRETERM_Thal_TC_norm,TERM_Thal_TC_norm);
% PREvsTERM_diff_thal = stats.tstat;

figure('Position',[1 49 1382 747])
for i = 1:6
    InDATA =  DATA2(:,i);
    subplot(2,3,i)
    s = scatterfit(InDATA,PREvsTERM_diff_thal,40,lines(1),[],2);
    xlabel(DATA2_LABEL{i})
ylabel({"Thalamic preterm-term","{\itF} statistic"})
    set(gca, 'FontSize', 16)
s.MarkerEdgeColor = [0 0 0];
end

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_thal_diff_scatts.png'],'-dpng','-r300')
