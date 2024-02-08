SurfData = zeros(length(medwallmask),1);

TERMvsPRE_Cort_norm_diff = mean(PRETERM_Cort_TC_norm{3})-mean(TERM_Cort_TC_norm{3});

SurfData(medwallmask) = PREvsTERM_diff;

plotRange = max(abs(SurfData));

SurfData(SurfData==0)=NaN;

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),'Mean difference',[-plotRange plotRange]);

ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),TYPENAME{TYPE},[-6 6]);

TERMvsPRE_Thal_norm_diff = mean(PRETERM_Thal_TC_norm)-mean(TERM_Thal_TC_norm);

up = max(abs(TERMvsPRE_Thal_norm_diff));
PlotThalGradientSlices(TERMvsPRE_Thal_norm_diff,vox_coords,cmocean('Balance'),['Preterm relative to term difference'],4,[-up up]);


% gL = gifti('week-40_hemi-right_space-dhcpSym_dens-32k_sphere.surf.gii');
% vertices = double(gL.vertices);
% spins = SpinVerts(vertices,100);