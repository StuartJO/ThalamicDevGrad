

TERMavgNorm = BF_NormalizeMatrix(TERMavg,'scaledSigmoid');

PRETERMavgNorm = BF_NormalizeMatrix(PRETERMavg,'scaledSigmoid');

scatter(TERMavgNorm(:),PRETERMavgNorm(:))

pvt = mean(PRETERMavg)-mean(TERMavg);

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = pvt;
up = max(abs(pvt));
RdBlu_cmap= cmocean('Balance',1024);
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),'Diff',[-up up]);


pvt = mean(PRETERMavgNorm)-mean(TERMavgNorm);

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = pvt;
up = max(abs(pvt));
RdBlu_cmap= cmocean('Balance',1024);
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),'Diff',[-up up]);



scatter(,,36,PvT_cortDiff,'filled')


NegThalDiff = (PCin{1}.PvT_thalDiff.*PCin{1}.fdr_thal)<0;
PosThalDiff = (PCin{1}.PvT_thalDiff.*PCin{1}.fdr_thal)>0;

NegTERM = nanmean(TERMavgNorm(NegThalDiff,:));
NegPRETERM = nanmean(PRETERMavgNorm(NegThalDiff,:));

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = NegTERM;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,parula(1024),'Connectivity of term (neg)'); 

SurfData(medwallmask) = NegPRETERM;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,parula(1024),'Connectivity of preterm (neg)'); 


SurfData(medwallmask) = NegPRETERM-NegTERM;
limits_abs = max(abs(NegPRETERM-NegTERM));
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance',1024),'Connectivity of preterm-term (neg)',[-limits_abs limits_abs]); 

%%

PosTERM = mean(TERMavgNorm(PosThalDiff,:));
PosPRETERM = mean(PRETERMavgNorm(PosThalDiff,:));

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = PosTERM;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,parula(1024),'Connectivity of term (pos)'); 

SurfData(medwallmask) = PosPRETERM;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,parula(1024),'Connectivity of preterm (pos)'); 


SurfData(medwallmask) = PosPRETERM-PosTERM;
limits_abs = max(abs(PosPRETERM-PosTERM));
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance',1024),'Connectivity of preterm-term (pos)',[-limits_abs limits_abs]); 



SurfData(medwallmask) = mean(PRETERM_Cort_TC_PC_norm{1,4})-mean(TERM_Cort_TC_PC_norm{1,4});
limits_abs = max(abs(mean(PRETERM_Cort_TC_PC_norm{1,4})-mean(TERM_Cort_TC_PC_norm{1,4})));
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance',1024),'Connectivity of preterm-term (pos)',[-limits_abs limits_abs]); 
