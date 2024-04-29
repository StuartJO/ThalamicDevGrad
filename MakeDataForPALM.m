load('./outputs/TermVsPretermTC_WEI_Thr100.mat')
load('UsedSubData.mat')

medmaskgii = gifti('week-40_hemi-left_space-dhcpSym_dens-32k_desc-medialwallsymm_mask.shape.gii');

medwallmask = logical(medmaskgii.cdata);

DATAMAT = zeros(102,length(medwallmask));
DATAMAT(1:51,medwallmask)=TERM_Cort_TC_norm;

DATAMAT(52:102,medwallmask)=PRETERM_Cort_TC_norm;

DESIGNMAT = zeros(102,4);

MASKMAT = double(medwallmask');

DESIGNMAT(1:51,1) = 1;
DESIGNMAT(52:102,2) = 1;

DESIGNMAT(:,3) = detrend([thal_sub_meta.scan_age(MatchedTERM); Preterm_scan_age],'constant');
DESIGNMAT(:,4) = detrend([SexInd(MatchedTERM); Preterm_sex],'constant');

mkdir ./PALM_results

writematrix(DESIGNMAT,'./PALM_results/design.csv')
writematrix(DATAMAT,'./outputs/data.csv')
writematrix(MASKMAT,'./PALM_results/mask.csv')

%CONMAT=[1 -1 0 0;-1 1 0 0;0 0 1 0;0 0 -1 0;0 0 0 1;0 0 0 -1];
CONMAT=[1 -1 0 0];
writematrix(CONMAT,'./PALM_results/contrasts.csv')

%% Run the following in bash

% awk 'BEGIN { FS=","; OFS=" " } { $1=$1; print $0 }' design.csv > design.txt
% awk 'BEGIN { FS=","; OFS=" " } { $1=$1; print $0 }' contrasts.csv > contrasts.txt
% module load fsl #May not need this depending on system configuration
% Text2Vest design.txt design.mat
% Text2Vest contrasts.txt design.con


% For deets see: https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM/CreatingDesignMatricesByHand

%% Once the above is made, run the following in MATLAB

% addpath ./PALM_results/
% addpath ./PALM-master/
%palm -i data.csv -d design.mat -t design.con -o ./PALM_results/TermVsPreterm -T -tfce2D -fdr -twotail -s week-40_hemi-left_space-dhcpSym_dens-32k_wm.surf.gii -logp 
%palm -i data.csv -d design.mat -t design.con -o ./PALM_results/TermVsPreterm_NO_TFCE -fdr -twotail -s week-40_hemi-left_space-dhcpSym_dens-32k_wm.surf.gii -logp
