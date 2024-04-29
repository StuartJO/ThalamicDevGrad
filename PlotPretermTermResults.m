Thr = 100;
Weighted=1;

if Weighted == 0
WEITYPE = 'Unweighted';
WEITYPE_ABBREV='_UNWEI';
elseif Weighted == 1
    WEITYPE = 'Weighted';
    WEITYPE_ABBREV='_WEI';
end

PCAOUTPUTFIG_DIR = ['./figures/Thr',num2str(Thr),'/',WEITYPE];

load(['.\outputs\TermVsPretermTC',WEITYPE_ABBREV,'_Thr',num2str(Thr),'.mat'])

% Get the TFCE result
T = dlmread('TermVsPreterm_tfce_tstat.csv');
% Get the TFCE p values
logp = dlmread('TermVsPreterm_tfce_tstat_fwep.csv');
% Get the normal t-statistic
t = dlmread('TermVsPreterm_NO_TFCE_dat_tstat.csv');

% Read in the gifti data
medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');
medwallmask = logical(medmaskgii.cdata);
gL = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = gL.vertices;
surface.faces = gL.faces;

% Transform the p values from log10 to non-log

p = 10.^(-logp);

% Flip the t statistic so it makes more sense to me :)
SurfData=-1*t;

SurfData_pmasked = SurfData;

%SurfData(medwallmask)=zscore(SurfData(medwallmask));
SurfDataRange = [-max(abs(SurfData)) max(abs(SurfData))];
SurfData_pmasked(p>.05) = NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData_pmasked,cmocean('Balance'),"Preterm-Term {\itt}-statistic",SurfDataRange,'none');

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PretermVsTerm_PALM_result.png'],'-dpng','-r300')

SurfDataMasked = SurfData(medwallmask);

PRETERM_Pos_mean = mean(PRETERM_Cort_TC_norm(:,SurfDataMasked>0),2);
TERM_Pos_mean = mean(TERM_Cort_TC_norm(:,SurfDataMasked>0),2);
PRETERM_Neg_mean = mean(PRETERM_Cort_TC_norm(:,SurfDataMasked<0),2);
TERM_Neg_mean = mean(TERM_Cort_TC_norm(:,SurfDataMasked<0),2);

lines_cmap = lines(4);
figure
% swarmdataX = [ones(1,51) ones(1,51).*2 ones(1,51).*4 ones(1,51).*5];
% swarmdataY = [TERM_Neg_mean; PRETERM_Neg_mean; TERM_Pos_mean; PRETERM_Pos_mean];
% swarmdataC = [ones(1,51) ones(1,51).*2 ones(1,51) ones(1,51).*2];
% swarmchart(swarmdataX,swarmdataY,50,swarmdataC,'filled')
% colormap(lines(2))
swarmdataX= [ones(1,51) ones(1,51).*4];
swarmdataY = [TERM_Neg_mean; TERM_Pos_mean];
s = swarmchart(swarmdataX,swarmdataY,50,lines_cmap(3,:),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s.XJitterWidth = 1;
hold on
swarmdataX= [ones(1,51).*2 ones(1,51).*5];
swarmdataY = [PRETERM_Neg_mean; PRETERM_Pos_mean];
s1 = swarmchart(swarmdataX,swarmdataY,50,lines_cmap(4,:),'filled','MarkerFaceAlpha',0.5,'MarkerEdgeAlpha',0.5);
s1.XJitterWidth = 1;
xticks([1.5 4.5])
xticklabels({'Term>Preterm','Term<Preterm'})



% Get the jitted position. Useful for me!
raw_s = struct(s);
raw_data = raw_s.XYZJittered;

jit = raw_data(1:51,1);
jit1 = raw_data(52:102,1);
meanline_xpos = [min(jit) min(jit1); max(jit) max(jit1);];
meanline_ypos = [mean(TERM_Neg_mean) mean(TERM_Pos_mean);mean(TERM_Neg_mean) mean(TERM_Pos_mean)];
meanplotTerm = plot(meanline_xpos,meanline_ypos,'LineWidth',4,'Color',lines_cmap(3,:));

raw_s1 = struct(s1);
raw_data = raw_s1.XYZJittered;

jit = raw_data(1:51,1);
jit1 = raw_data(52:102,1);
meanline_xpos = [min(jit) min(jit1); max(jit) max(jit1);];
meanline_ypos = [mean(PRETERM_Neg_mean) mean(PRETERM_Pos_mean);mean(PRETERM_Neg_mean) mean(PRETERM_Pos_mean)];
meanplotPreterm = plot(meanline_xpos,meanline_ypos,'LineWidth',4,'Color',lines_cmap(4,:));
lgd = legend([meanplotTerm(1) meanplotPreterm(1)],'Term','Preterm','Orientation','horizontal','Location','northoutside');
ylabel('Mean cluster connectivity')
set(gca,'FontSize',20)

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PretermVsTerm_PALM_cluster_result.png'],'-dpng','-r300')


%%
load('.\outputs\SpunVerts.mat')

PC1 = load(['./outputs/VERT',WEITYPE_ABBREV,'_PC1_Thr',num2str(Thr),'_Results.mat']);
PC2 = load(['./outputs/VERT',WEITYPE_ABBREV,'_PC2_Thr',num2str(Thr),'_Results.mat']);
PC3 = load(['./outputs/VERT',WEITYPE_ABBREV,'_PC3_Thr',num2str(Thr),'_Results.mat']);

DATA(:,1) = TemplateCoeffz(:,1);
DATA(:,2) = TemplateCoeffz(:,2);
DATA(:,3) = TemplateCoeffz(:,3);

DATA(:,4) = PC1.PC_cort_age_corr;
DATA(:,5) = PC2.PC_cort_age_corr;
DATA(:,6) = PC3.PC_cort_age_corr;

DATA_LABEL = {'PC1 template loading','PC2 template loading','PC3 template loading','PC1 age loading corr.','PC2 age loading corr.','PC3 age loading corr.'};

figure('Position',[1 49 1382 747])

for i = 1:6

InDATA = NaN(length(medwallmask),1);    
InDATA(medwallmask) = DATA(:,i);

p_perm = perm_sphere_p(InDATA,SurfData',spins','Pearson');

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
ylabel("Preterm-Term {\itt}-statistic")

set(gca, 'FontSize', 16)
s.MarkerEdgeColor = [0 0 0];

end

print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'pretermVSterm_cort_diff_scatts.png'],'-dpng','-r300')

%%

PretermIND = ones(length(PretermTermScanInd),1);
TermIND = ones(length(MatchedTermInd),1)*2;
GrpInd = [PretermIND;TermIND];

AgeAtScan = [PretermTermScan_scan_age ;MatchedTerm_scan_age];
NeonateSex = [PretermTermScan_sex;MatchedTerm_sex];

TERMvsPRE_ANVOCA = zeros(3,6);

for i = 1:3

DATA_MAT(:,1) = GrpInd;
DATA_MAT(:,2) = [AllExpl(PretermTermScanInd,i);AllExpl(MatchedTermInd,i)];
DATA_MAT(:,3) = AgeAtScan;
DATA_MAT(:,4) = NeonateSex;

[p,tbl] = anovan([AllExpl(PretermTermScanInd,i);AllExpl(MatchedTermInd,i)],{GrpInd,AgeAtScan,NeonateSex},"Continuous",[2],"Varnames",["TermVsPreterm","Age at scan","sex"],'display','off');
TERMvsPRE_ANVOCA(i,2) = p(1);
TERMvsPRE_ANVOCA(i,1) = tbl{2,6};
TERMvsPRE_ANVOCA(i,[3 4]) = [mean(AllExpl(PretermTermScanInd,i)) mean(AllExpl(MatchedTermInd,i))];
TERMvsPRE_ANVOCA(i,[5 6]) = [std(AllExpl(PretermTermScanInd,i)) std(AllExpl(MatchedTermInd,i))];
end
