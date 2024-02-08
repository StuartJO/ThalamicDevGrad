

PretermIND = ones(length(PRETERM),1);
TermIND = ones(length(MatchedSUB),1)*2;
GrpInd = [PretermIND;TermIND];
AgeAtScan = [Preterm_scan_age ;thal_sub_meta.scan_age(MatchedSUB)];
ancova_thal = zeros(Nseed,1);
F_thal = zeros(Nseed,1);
for i = 1:Nseed
    [p,tbl,stats] = anovan([PC_thal_all(i,PRETERM)';PC_thal_all(i,MatchedSUB)'],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
ancova_thal(i,:) = p(1);
F_thal(i,:) = tbl{2,6};
end

fdr_thal = BF_FDR(ancova_thal,.05);
thal_diff = mean(PC_thal_all(:,PRETERM)')-mean(PC_thal_all(:,MatchedSUB)');
thal_diff_sign = sign(thal_diff)';

val = F_thal.*thal_diff_sign.*fdr_thal;


PC=2;
figure
scatter(mean(PC_thal_all(:,MatchedSUB)'),mean(PC_thal_all(:,PRETERM)'),36,val,'filled')
colormap(cmocean('Balance'))
up = max(abs(val));
caxis([-up up])
xlabel(['Term PC',num2str(PC),' scores'])
ylabel(['Preterm PC',num2str(PC),' scores'])
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' score difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)

PlotThalGradientSlices(val,vox_coords,cmocean('Balance'),['Preterm relative to term PC',num2str(PC),' score difference'],4,[-up up]);



[h1,p1,~,t1]= ttest2((PC_cort_all(:,PRETERM)'),(PC_cort_all(:,MatchedSUB)'));

fdr = BF_FDR(p1,.05)';


PretermIND = ones(length(PRETERM),1);
TermIND = ones(length(MatchedSUB),1)*2;
GrpInd = [PretermIND;TermIND];
AgeAtScan = [Preterm_scan_age ;thal_sub_meta.scan_age(MatchedSUB)];
ancova_cort = zeros(Nverts_nomed,1);
F = zeros(Nverts_nomed,1);
for i = 1:Nverts_nomed
    [p,tbl,stats] = anovan([PC_cort_all(i,PRETERM)';PC_cort_all(i,MatchedSUB)'],{GrpInd,AgeAtScan},"Continuous",2,"Varnames",["TermVsPreterm","Age at scan"],'display','off');
ancova_cort(i,:) = p(1);
F(i,:) = tbl{2,6};
end

fdr_cort = BF_FDR(ancova_cort,.05);

cort_diff = mean(PC_cort_all(:,PRETERM)')-mean(PC_cort_all(:,MatchedSUB)');

diff_sign = sign(cort_diff)';

figure

val = F.*diff_sign.*fdr_cort;

scatter(mean(PC_cort_all(:,MatchedSUB)'),mean(PC_cort_all(:,PRETERM)'),10,val,'filled')
colormap(cmocean('Balance'))
up = max(abs(val));
caxis([-up up])
xlabel(['Term PC',num2str(PC),' loadings'])
ylabel(['Preterm PC',num2str(PC),' loadings'])
c = colorbar;
c.Label.String = {['Preterm relative to term PC',num2str(PC)],' loadings difference'};
c.Label.FontSize = 14;
set(gca,'FontSize',18)


val = F.*diff_sign;
SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = val;
up = max(abs(val));
RdBlu_cmap= cmocean('Balance',1024);

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),'F',[-up up]); 

val(~fdr_cort) = NaN;
SurfData(medwallmask) = val;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),'F',[-up up]); 

figure
scatter(F.*diff_sign,mean(PRETERM_Cort_TC_norm)-mean(TERM_Cort_TC_norm))

