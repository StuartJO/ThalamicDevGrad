load('SeedsThr.mat','ThalConnMean')
SeedThr = ThalConnMean>=100;
vox_coords_full = dlmread('thal_seed_1.75mm_vox_coords.txt');
vox_coords = vox_coords_full(SeedThr,:);

TERMavgNorm = BF_NormalizeMatrix(TERMavg,'scaledSigmoid');

PRETERMavgNorm = BF_NormalizeMatrix(PRETERMavg,'scaledSigmoid');

sig = (sign(PvT_thalDiff).*fdr_thal)==1;
issig =  find(sig==1);
[~,ord] = sort(PvT_thalDiff,'descend');

%[~,ord] = sort(TrainScorez(:,1),'descend');

%ord(ismember(ord,issig))=[];

%InClust = ismember(ord,find(sig));



Nsig = sum(sig);

Nseed = length(sig);

Ntests = length(ord)-(Nsig-1);



% TermMean = mean(mean(TERMavgNorm(sig,:)));
% PretermMean = mean(mean(PRETERMavgNorm(sig,:)));
TermMean = mean(TERMavgNorm(sig,:),'all');
PretermMean = mean(PRETERMavgNorm(sig,:),'all');


notsig = find(sig==0);
Nnotsig = length(notsig);
Nperm = Ntests;
TermMean_=zeros(Nperm,1);
PretermMean_=zeros(Nperm,1);
PermD=zeros(Nperm,1);
meanD = mean(pdist(vox_coords(sig,:)));
clustOverlap = zeros(Nperm,1);
for i = 1:Nperm
Perm = ord(i:i+(Nsig-1));
PermD(i) = mean(pdist(vox_coords(Perm,:)));
TermMean_(i) = mean(TERMavgNorm(Perm,:),'all');
PretermMean_(i) = mean(PRETERMavgNorm(Perm,:),'all');
clustOverlap(i) = sum(ismember(Perm,issig))/Nsig;
end
histogram(PretermMean_-TermMean_)
ylimits = ylim;
hold on
plot([PretermMean-TermMean PretermMean-TermMean],ylimits,'r')
figure
plot(PretermMean_-TermMean_)
% plotyy(1:Nperm,PretermMean_-TermMean_,1:Nperm,clustOverlap)

% notsig = find(sig==0);
% Nnotsig = length(notsig);
% Nperm = 10000;
% TermMean_=zeros(Nperm,1);
% PretermMean_=zeros(Nperm,1);
% PermD=zeros(Nperm,1);
% meanD = mean(pdist(vox_coords(sig,:)));
% for i = 1:Nperm
% perm = randperm(Nnotsig,Nsig);
% Perm = notsig(perm);
% PermD(i) = mean(pdist(vox_coords(Perm,:)));
% TermMean_(i) = mean(TERMavgNorm(Perm,:),'all');
% PretermMean_(i) = mean(PRETERMavgNorm(Perm,:),'all');
% end
% histogram(PretermMean_-TermMean_)
% ylimits = ylim;
% hold on
% plot([PretermMean-TermMean PretermMean-TermMean],ylimits,'r')

% 
gL = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_vinflated.surf.gii');
surface.vertices = gL.vertices;
surface.faces = gL.faces;

SurfData = zeros(length(medwallmask),1);
SurfData(medwallmask) = TrainCoeffz(:,8);
up = max(abs(SurfData));
RdBlu_cmap= cmocean('Balance',1024);
ExampleSurfacePlotFunction(surface,double(medwallmask),SurfData,cmocean('Balance'),'Diff',[-up up]);

% 

[Coeff1,Score1,~,~,Explained1] = pca(TERMavgNorm);
[Coeff2,Score2,~,~,Explained2] = pca(PRETERMavgNorm);

Norm = BF_NormalizeMatrix(TrainDataAvg,'scaledSigmoid');
Norm(isnan(Norm)) = 0;
[TrainCoeff,TrainScore,~,~,TrainExplained] = pca(Norm);

[aligned_pvt, xfms_pvt] = procrustes_alignment({Score1,Score2},'reference',TrainScore);

Coeff1a = Coeff1*xfms_pvt{1};
Coeff2a = Coeff2*xfms_pvt{2};

T_PCs = zscore(aligned_pvt{1});
P_PCs = zscore(aligned_pvt{2});

for i = 1:10
    subplot(2,5,i)
pc = i;
di = P_PCs(:,pc)-T_PCs(:,pc);
scatter(T_PCs(:,pc),P_PCs(:,pc),20,di,'filled')
colormap(cmocean('Balance'))
limit = max(abs(di));
caxis([-limit limit])       
end