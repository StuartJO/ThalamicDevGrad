


Y = dlmread('TermVsPreterm_tfce_tstat_c1.csv');
X = dlmread('TermVsPreterm_tfce_tstat_c2.csv');
logp = dlmread('TermVsPreterm_tfce_tstat_fdrp_c2.csv');

p = 10.^(-logp);

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');
medwallmask = logical(medmaskgii.cdata);
gL = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = gL.vertices;
surface.faces = gL.faces;

% SurfData = zeros(length(medwallmask),1);
% SurfData(medwallmask) = TERMvsPRE_Cort_norm(:,1).*sign(TERMvsPRE_Cort_norm(:,3));
% plotRange = max(abs(SurfData));
% SurfData(SurfData==0)=NaN;

SurfData = X;
SurfData(p>.05)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"Preterm-term {\itF} statistic",[min(X) max(X)],'none');


Pre_Term_diff = mean(PRETERM_Cort_TC_norm)-mean(TERM_Cort_TC_norm);
Pre_Term_diff_sign = sign(Pre_Term_diff);
SurfData = X;
SurfData(medwallmask)=SurfData(medwallmask).*Pre_Term_diff_sign;
SurfDataRange = [min(SurfData) max(SurfData)];
SurfData(p>.05)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"Preterm-term {\itF} statistic",SurfDataRange,'none');

Pre_Term_diff = TERMvsPRE_Cort_norm(:,1);

scatter(X(medwallmask),Pre_Term_diff)

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"PALM result",[min(X) max(X)],'none');

SurfData(medwallmask)=Pre_Term_diff;
SurfDataRange = [min(Pre_Term_diff) max(Pre_Term_diff)];
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"Preterm-term {\itF} statistic",SurfDataRange,'none');
%%
T_P = dlmread('TermVsPreterm_tfce_tstat_c1.csv');
logp = dlmread('TermVsPreterm_tfce_tstat_fwep_c1.csv');
p = 10.^(-logp);
Pre_Term_diff = mean(PRETERM_Cort_TC_norm)-mean(TERM_Cort_TC_norm);
Pre_Term_diff_sign = sign(Pre_Term_diff);
SurfData = T_P;
SurfData(medwallmask)=zscore(SurfData(medwallmask).*Pre_Term_diff_sign);
SurfDataRange = [-max(abs(SurfData)) max(abs(SurfData))];
SurfData(p>.05)=NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"Preterm-term TFCE adjusted {\itt} statistic",SurfDataRange,'none');

print('PretermTerm_PALM_fwer_surface_plot.png','-dpng','-r300')

SurfDataMasked = SurfData(medwallmask);

PRETERM_Pos_mean = mean(mean(PRETERM_Cort_TC_norm(:,SurfDataMasked>0)));
TERM_Pos_mean = mean(mean(TERM_Cort_TC_norm(:,SurfDataMasked>0)));
PRETERM_Neg_mean = mean(mean(PRETERM_Cort_TC_norm(:,SurfDataMasked<0)));
TERM_Neg_mean = mean(mean(TERM_Cort_TC_norm(:,SurfDataMasked<0)));
bar([TERM_Neg_mean PRETERM_Neg_mean ; TERM_Pos_mean PRETERM_Pos_mean])
legend('Term','Preterm')

PRETERM_Pos_mean = mean(PRETERM_Cort_TC_norm(:,SurfDataMasked>0),2);
TERM_Pos_mean = mean(TERM_Cort_TC_norm(:,SurfDataMasked>0),2);
PRETERM_Neg_mean = mean(PRETERM_Cort_TC_norm(:,SurfDataMasked<0),2);
TERM_Neg_mean = mean(TERM_Cort_TC_norm(:,SurfDataMasked<0),2);

figure
% swarmdataX = [ones(1,51) ones(1,51).*2 ones(1,51).*4 ones(1,51).*5];
% swarmdataY = [TERM_Neg_mean; PRETERM_Neg_mean; TERM_Pos_mean; PRETERM_Pos_mean];
% swarmdataC = [ones(1,51) ones(1,51).*2 ones(1,51) ones(1,51).*2];
% swarmchart(swarmdataX,swarmdataY,50,swarmdataC,'filled')
% colormap(lines(2))
swarmdataX= [ones(1,51) ones(1,51).*4];
swarmdataY = [TERM_Neg_mean; TERM_Pos_mean];
s = swarmchart(swarmdataX,swarmdataY,50,'filled');
s.XJitterWidth = 1;
hold on
swarmdataX= [ones(1,51).*2 ones(1,51).*5];
swarmdataY = [PRETERM_Neg_mean; PRETERM_Pos_mean];
s1 = swarmchart(swarmdataX,swarmdataY,50,'filled');
s1.XJitterWidth = 1;
legend('Term','Preterm')
xticks([1.5 4.5])
xticklabels({'Term>Preterm','Term<Preterm'})
% swarmchart(ones(1,51),TERM_Neg_mean,50,'filled')
% hold on
% swarmchart(ones(1,51).*2,PRETERM_Neg_mean,50,'filled')
% swarmchart(ones(1,51).*4,TERM_Pos_mean,50,'filled','F')
% swarmchart(ones(1,51).*5,PRETERM_Pos_mean,50,'filled')
% 
% hold off
%%

TgtP = dlmread('TermVsPreterm_OneTail_detrend_tfce_tstat_c1.csv');
TltP = dlmread('TermVsPreterm_OneTail_detrend_tfce_tstat_c2.csv');

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');
medwallmask = logical(medmaskgii.cdata);
gL = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = gL.vertices;
surface.faces = gL.faces;

SurfData=TltP-TgtP;
SurfData(medwallmask)=zscore(SurfData(medwallmask));
SurfDataRange = [-max(abs(SurfData)) max(abs(SurfData))];
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"PALM result",SurfDataRange,'none');

[h,p,ci,stats] = ttest2(PRETERM_Cort_TC_norm,TERM_Cort_TC_norm);

T = stats.tstat;

SurfData=TltP-TgtP;

scatter(T,SurfData(medwallmask))

%%

TgtP = dlmread('TermVsPreterm_NO_TFCE_dat_tstat_c1.csv');
TltP = dlmread('TermVsPreterm_NO_TFCE_dat_tstat_c2.csv');