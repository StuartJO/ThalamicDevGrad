
makefig = 1;
saveoutput=1;

Weighted=1;

if Weighted == 0
WEITYPE = 'Unweighted';
WEITYPE_ABBREV='_UNWEI';
elseif Weighted == 1
    WEITYPE = 'Weighted';
    WEITYPE_ABBREV='_WEI';
end

PCAOUTPUTFIG_DIR = ['./figures/Thr',num2str(Thr),'/',WEITYPE];

mkdir(PCAOUTPUTFIG_DIR)

for PC=1:3

for i = 1:length(SUB)

PC_thal_all(:,i) = zscore(aligned{i}(:,PC));
PC_cort_all(:,i) =  zscore((coeff_aligned_all{i}(:,PC)));

end

PC_thal = PC_thal_all(:,NonTemplateInd);
PC_cort = PC_cort_all(:,NonTemplateInd);

TemplateCoeffz_cort = zeros(Nverts,Nseed-1);
TemplateCoeffz_cort(medwallmask,:) = TemplateCoeffz;

gL = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = gL.vertices;
surface.faces = gL.faces;

CorrType = 'Pearson';

[PC_corr_wAvg,PC_corr_wAvg_p] = corr(TemplateCoeffz(:,PC),PC_cort,'Type',CorrType);

[PC_corr_wAvg_thal,PC_corr_wAvg_thal_p] = corr(TemplateScorez(:,PC),PC_thal,'Type',CorrType);

[PC_age_corr,PC_age_corr_p] = corr(NonTemplate_scan_age,PC_cort','Type',CorrType);

CorrRange = max(abs(PC_age_corr(:)));
up = ceil(CorrRange * 10) / 10;

fdr = BF_FDR(PC_age_corr_p,.05)';

PC_age_corr_cort = zeros(Nverts,1);
PC_age_corr_cort(medwallmask) = PC_age_corr.*fdr;
%PC_age_corr_cort(medwallmask) = PC_age_corr; 
PC_age_corr_cort(PC_age_corr_cort==0)=NaN;

ExampleSurfacePlotFunction(surface,double(medwallmask),PC_age_corr_cort,cmocean('Balance'),['PC',num2str(PC),' loading correlation with scan age'],[-up up],'none');    

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_loading_corr_scan_age_region.png'],'-dpng','-r300')
end

ExampleSurfacePlotFunction(surface,double(medwallmask),TemplateCoeffz_cort(:,PC),turbo(256),['Term template PC',num2str(PC),' loading'],[min(TemplateCoeffz_cort(:,PC)) max(TemplateCoeffz_cort(:,PC))],'none');    

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_Group_PC',num2str(PC),'_loading.png'],'-dpng','-r300')
end

figure %figure('Position',[189   318   946   688])
s = scatterfit(TemplateCoeffz(:,PC),PC_age_corr,36,PC_age_corr.*fdr,[],2);
s.MarkerFaceAlpha=1;
colormap(cmocean('Balance'))
maxCorr = max(abs(PC_age_corr.*fdr));
caxis([-maxCorr maxCorr])
xlabel(['Term template PC',num2str(PC),' loading'])
ylabel({['PC',num2str(PC),' loading correlation'],'with scan age'})
set(gca,'FontSize',18)
c = colorbar;
c.Label.String = {['PC',num2str(PC),' loading correlation with scan age']};
c.Label.FontSize = 14;

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_Group_PC',num2str(PC),'_loading_vs_scan_age_corr.png'],'-dpng','-r300')
end

[PC_thal_age_corr,PC_thal_age_corr_p] = corr(NonTemplate_scan_age,PC_thal','Type',CorrType);
fdr_thal = BF_FDR(PC_thal_age_corr_p,.05)';
CorrRange = max(abs(PC_thal_age_corr(:)));
up = ceil(CorrRange * 10) / 10;

vox_coords_full = dlmread('thal_seed_1.75mm_vox_coords.txt');

vox_coords = vox_coords_full(SeedThr,:);

PlotThalGradientSlices((TemplateScorez(:,PC)),vox_coords,turbo(256),['Term template PC',num2str(PC),' score'],4);
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_Group_PC',num2str(PC),'_score.png'],'-dpng','-r300')
end

ThalData = PC_thal_age_corr.*fdr_thal;
ThalData(~fdr_thal)=NaN;
PlotThalGradientSlices(ThalData,vox_coords,cmocean('Balance'),['PC',num2str(PC),' score correlation with scan age'],4,[-up up]);
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_score_corr_scan_age_region.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(TemplateScorez(:,PC),PC_thal_age_corr,36,PC_thal_age_corr.*fdr_thal,[],2);
s.MarkerFaceAlpha=1;
colormap(cmocean('Balance'))
caxis([-up up])
xlabel(['Term template PC',num2str(PC),' score'])
ylabel({['PC',num2str(PC),' score correlation'],'with scan age'})
set(gca,'FontSize',18)
c = colorbar;
c.Label.String = {['PC',num2str(PC),' score correlation with scan age']};
c.Label.FontSize = 14;
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_Group_PC',num2str(PC),'_score_vs_scan_age_corr.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(NonTemplate_scan_age,PC_corr_wAvg);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({['Correlation with PC',num2str(PC)],'term template (loading)'})
set(gca,'FontSize',18)
if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_loading_corr_scan_age.png'],'-dpng','-r300')
end
figure %figure('Position',[189   318   946   688])
s = scatterfit(NonTemplate_scan_age,PC_corr_wAvg_thal);
s.MarkerFaceAlpha=1;
xlabel('Scan age (weeks)')
ylabel({['Correlation with PC',num2str(PC)],'term template (score)'})
set(gca,'FontSize',18)

s.MarkerFaceAlpha=1;

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_score_corr_scan_ageALPHA.png'],'-dpng','-r300')
end

figure('Position',[62   300   1300   335])
for i = 1:3
    subplot(1,3,i)
[s,t] = scatterfit(vox_coords(:,i),PC_thal_age_corr,36,TemplateScorez(:,PC));
s.MarkerFaceAlpha=1;
cardDir={'Medial-Lateral','Anterior-posterior','Ventral-Dorsal'};
xlabel(cardDir{i})
ylabel({['PC',num2str(PC),' score correlation'],'with scan age'})
set(gca,'FontSize',12)
t.FontSize = 12;
tString = t.String;
t.String = {[tString{1},' ',tString{2}]};
t.Position(1) = sum(xlim)/2;
t.HorizontalAlignment='center';
ylimits = ylim;
t.Position(2) = find_point_on_line(ylimits(1),ylimits(2),.95);
colormap(turbo(256))
end
c = colorbar;
c.Position = [0.927282051527806	0.146778037162555	0.0151076923076924	0.778042971020951];
c.Label.String = ['Term template PC',num2str(PC),' loading'];

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_score_CardDirCorr.png'],'-dpng','-r300')
end

figure('Position',[62   300   1300   335])

for i = 1:3
    subplot(1,3,i)
[s,t] = scatterfit(vox_coords(:,i),TemplateScorez(:,PC),36,lines(1),[],2);
s.MarkerFaceAlpha=.5;
cardDir={'Medial-Lateral','Anterior-posterior','Ventral-Dorsal'};
xlabel(cardDir{i})
ylabel(['PC',num2str(PC),' score'])
set(gca,'FontSize',12)
% t.FontSize = 12;
% tString = t.String;
% t.String = {[tString{1},' ',tString{2}]};
% t.Position(1) = sum(xlim)/2;
% t.HorizontalAlignment='center';
% ylimits = ylim;
% t.Position(2) = find_point_on_line(ylimits(1),ylimits(2),.95);
% colormap(turbo(256))
end
% c = colorbar;
% c.Position = [0.927282051527806	0.146778037162555	0.0151076923076924	0.778042971020951];
% c.Label.String = ['Term template PC',num2str(PC),' loading'];

if makefig
print([PCAOUTPUTFIG_DIR,'/VERT',WEITYPE_ABBREV,'_PC',num2str(PC),'_score_CardDir.png'],'-dpng','-r300')
end

close all

end