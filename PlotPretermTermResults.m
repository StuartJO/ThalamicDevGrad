T = dlmread('TermVsPreterm_TFCE_tfce_tstat.csv');
P = dlmread('TermVsPreterm_TFCE_tfce_tstat_fewp.csv');
t = dlmread('TermVsPreterm_NO_TFCE_dat_tstat.csv');
medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');
medwallmask = logical(medmaskgii.cdata);
gL = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = gL.vertices;
surface.faces = gL.faces;

SurfData=t;
%SurfData(medwallmask)=zscore(SurfData(medwallmask));
SurfDataRange = [-max(abs(SurfData)) max(abs(SurfData))];
SurfData(P<.05) = NaN;
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),"PALM result",SurfDataRange,'none');
