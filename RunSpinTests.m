
g = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_sphere.surf.gii');
vertices = g.vertices;
spins = SpinVerts(vertices, 1000);
save('./outputs/SpunVerts.mat','spins')