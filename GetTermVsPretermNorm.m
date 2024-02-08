
if Weighted == 0
    suf = '_wei';
elseif Weighted == 1
    suf = '';
end

for PC = 1:2
PCin{PC} = load(['VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_Thr',num2str(Thr),'_Results.mat']);
end


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

for PC = 1:2

for i = 1:4

    PRETERM_Cort_TC_PC{PC,i} = zeros(length(PRETERM),nCortVerts);
    PRETERM_Cort_TC_PC_norm{PC,i} = zeros(length(PRETERM),nCortVerts);
    TERM_Cort_TC_PC{PC,i} = zeros(length(MatchedTERM),nCortVerts);
    TERM_Cort_TC_PC_norm{PC,i} = zeros(length(MatchedTERM),nCortVerts);

end

end

%PRETERM_SIGTC = zeros(length(SigSeeds),28766,length(PRETERM));

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

for PC = 1:2

SigSeeds = find(PCin{PC}.fdr_thal)';
NonSigSeeds = find(~PCin{PC}.fdr_thal)';
NegThalDiff = (PCin{PC}.PvT_thalDiff.*PCin{PC}.fdr_thal)<0;
PosThalDiff = (PCin{PC}.PvT_thalDiff.*PCin{PC}.fdr_thal)>0;


PRETERM_Cort_TC_PC{PC,1}(i,:) = sum(tc(SigSeeds,:));
PRETERM_Cort_TC_PC{PC,2}(i,:) = sum(tc(NonSigSeeds,:));

PRETERM_Cort_TC_PC{PC,3}(i,:) = sum(tc(NegThalDiff,:));
PRETERM_Cort_TC_PC{PC,4}(i,:) = sum(tc(PosThalDiff,:));

PRETERM_Cort_TC_PC_norm{PC,1}(i,:) = sum(norm(SigSeeds,:));
PRETERM_Cort_TC_PC_norm{PC,2}(i,:) = sum(norm(NonSigSeeds,:));

PRETERM_Cort_TC_PC_norm{PC,3}(i,:) = sum(norm(NegThalDiff,:));
PRETERM_Cort_TC_PC_norm{PC,4}(i,:) = sum(norm(PosThalDiff,:));

end

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

for PC = 1:2

SigSeeds = find(PCin{PC}.fdr_thal)';
NonSigSeeds = find(~PCin{PC}.fdr_thal)';
NegThalDiff = (PCin{PC}.PvT_thalDiff.*PCin{PC}.fdr_thal)<0;
PosThalDiff = (PCin{PC}.PvT_thalDiff.*PCin{PC}.fdr_thal)>0;

TERM_Cort_TC_PC{PC,1}(i,:) = sum(tc(SigSeeds,:));
TERM_Cort_TC_PC{PC,2}(i,:) = sum(tc(NonSigSeeds,:));

TERM_Cort_TC_PC{PC,3}(i,:) = sum(tc(NegThalDiff,:));
TERM_Cort_TC_PC{PC,4}(i,:) = sum(tc(PosThalDiff,:));

TERM_Cort_TC_PC_norm{PC,1}(i,:) = sum(norm(SigSeeds,:));
TERM_Cort_TC_PC_norm{PC,2}(i,:) = sum(norm(NonSigSeeds,:));

TERM_Cort_TC_PC_norm{PC,3}(i,:) = sum(norm(NegThalDiff,:));
TERM_Cort_TC_PC_norm{PC,4}(i,:) = sum(norm(PosThalDiff,:));

end

TERM_Thal_TC(i,:) = sum(tc,2);

TERM_Thal_TC_norm(i,:) = sum(norm,2);
end

save(['TermVsPretermTC',WEITYPE_ABBREV,'_Thr',num2str(Thr),'.mat'],'TERM_Cort_TC_PC_norm','PRETERM_Cort_TC_PC_norm','TERM_Cort_TC_PC','PRETERM_Cort_TC_PC','TERM_Thal_TC_norm','PRETERM_Thal_TC_norm','TERM_Thal_TC','PRETERM_Thal_TC','TERM_Cort_TC_norm','TERM_Cort_TC','PRETERM_Cort_TC_norm','PRETERM_Cort_TC')

TERMavg = TERMsum./length(MatchedTERM);
PRETERMavg = PRETERMsum./length(MatchedTERM);
clear TERMsum PRETERMsum
save('TermPretermAvg.mat','TERMavg','PRETERMavg')

nCortVerts = sum(medwallmask);
for PC = 2:-1:1

PretermIND = ones(length(PRETERM),1);
TermIND = ones(length(MatchedTERM),1)*2;
GrpInd = [PretermIND;TermIND];
AgeAtScan = [Preterm_scan_age ;thal_sub_meta.scan_age(MatchedTERM)];

TERMvsPRE_Cort_norm = cell(1,5);

for j = 1:4
    TERMvsPRE_Cort_norm{j} = zeros(nCortVerts,3);
        preterm_data = PRETERM_Cort_TC_PC_norm{PC,j};
        term_data =TERM_Cort_TC_PC_norm{PC,j};
for i = 1:nCortVerts
    [p,tbl] = anovan([preterm_data(:,i);term_data(:,i)],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
    TERMvsPRE_Cort_norm{j}(i,2) = p(1);
    TERMvsPRE_Cort_norm{j}(i,1) = tbl{2,6};
end
TERMvsPRE_Cort_norm{j}(:,3) = mean(preterm_data)-mean(term_data);
end
preterm_data = PRETERM_Cort_TC_norm;
term_data =TERM_Cort_TC_norm;
for i = 1:nCortVerts
    [p,tbl] = anovan([preterm_data(:,i);term_data(:,i)],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
    TERMvsPRE_Cort_norm{5}(i,2) = p(1);
    TERMvsPRE_Cort_norm{5}(i,1) = tbl{2,6};
end
TERMvsPRE_Cort_norm{5}(:,3) = mean(preterm_data)-mean(term_data);



TERMvsPRE_Cort = cell(1,5);
for j = 1:4
    TERMvsPRE_Cort{j} = zeros(nCortVerts,3);
        preterm_data = PRETERM_Cort_TC_PC{PC,j};
        term_data =TERM_Cort_TC_PC{PC,j};
for i = 1:nCortVerts
    [p,tbl] = anovan([preterm_data(:,i);term_data(:,i)],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
    TERMvsPRE_Cort{j}(i,2) = p(1);
    TERMvsPRE_Cort{j}(i,1) = tbl{2,6};
end
TERMvsPRE_Cort{j}(:,3) = mean(preterm_data)-mean(term_data);
end
preterm_data = PRETERM_Cort_TC;
term_data =TERM_Cort_TC;
for i = 1:nCortVerts
    [p,tbl] = anovan([preterm_data(:,i);term_data(:,i)],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
    TERMvsPRE_Cort{5}(i,2) = p(1);
    TERMvsPRE_Cort{5}(i,1) = tbl{2,6};
end
TERMvsPRE_Cort{5}(:,3) = mean(preterm_data)-mean(term_data);

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

TYPENAME = {'Sig. seeds difference (P-T)','Non-sig. seeds difference (P-T)','Negsig. seeds difference (P-T)','Possig. seeds difference (P-T)','Preterm vs term difference',};

for TYPE = 1:5

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort_norm{TYPE}(:,3).*(TERMvsPRE_Cort_norm{TYPE}(:,2)<.05);

plotRange = max(abs(SurfData));

SurfData(SurfData==0)=NaN;

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),TYPENAME{TYPE},[-plotRange plotRange]);

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_pretermVSterm_cort_',num2str(TYPE),'.png'],'-dpng','-r300')
end

end

thal_sig = TERMvsPRE_Thal(:,2)<.05;

up = max(abs(TERMvsPRE_Thal(:,3).*thal_sig));
PlotThalGradientSlices(TERMvsPRE_Thal(:,3).*thal_sig,vox_coords,cmocean('Balance'),['Preterm relative to term difference'],4,[-up up]);

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_pretermVSterm_thal_difference_corr.png'],'-dpng','-r300')
end


thal_sig = TERMvsPRE_Thal_norm(:,2)<1;

up = max(abs(TERMvsPRE_Thal_norm(:,3).*thal_sig));
PlotThalGradientSlices(TERMvsPRE_Thal_norm(:,3).*thal_sig,vox_coords,cmocean('Balance'),['Preterm relative to term difference'],4,[-up up]);

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_pretermVSterm_thal_difference_uncor.png'],'-dpng','-r300')
end

end

close all


SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort_norm{5}(:,3);
plotRange = max(abs(SurfData));
SurfData(SurfData==0)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),['Preterm relative to term difference'],[-plotRange plotRange]);
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_pretermVSterm_thal_difference_uncorr.png'],'-dpng','-r300')
end

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort_norm{5}(:,3).*(TERMvsPRE_Cort_norm{5}(:,2)<.05);
plotRange = max(abs(SurfData));
SurfData(SurfData==0)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),['Preterm relative to term difference'],[-plotRange plotRange]);
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_pretermVSterm_thal_difference_corr.png'],'-dpng','-r300')
end


SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TERMvsPRE_Cort{5}(:,1);
plotRange = max(abs(SurfData));
SurfData(SurfData==0)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),['Preterm relative to term difference'],[-plotRange plotRange]);



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

