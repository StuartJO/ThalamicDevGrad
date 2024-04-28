

DATAMAT = zeros(102,length(medwallmask));
DATAMAT(1:51,medwallmask)=TERM_Cort_TC_norm;

DATAMAT(52:102,medwallmask)=PRETERM_Cort_TC_norm;

DATAMAT_detrend = detrend(DATAMAT,'constant');

writematrix(DATAMAT_detrend,'data_detrend.csv')

DESIGNMAT = zeros(102,4);

MASKMAT = double(medwallmask');

DESIGNMAT(1:51,1) = 1;
DESIGNMAT(52:102,2) = 1;

DESIGNMAT(:,3) = detrend([thal_sub_meta.scan_age(MatchedTERM); Preterm_scan_age],'constant');
DESIGNMAT(:,4) = detrend([SexInd(MatchedTERM); Preterm_sex],'constant');

writematrix(DESIGNMAT,'design.csv')
writematrix(DATAMAT,'data.csv')
writematrix(MASKMAT,'mask.csv')

%CONMAT=[1 -1 0 0;-1 1 0 0;0 0 1 0;0 0 -1 0;0 0 0 1;0 0 0 -1];
CONMAT=[1 -1 0 0];
writematrix(CONMAT,'contrasts.csv')

%% Run the following in bash

% awk 'BEGIN { FS=","; OFS=" " } { $1=$1; print $0 }' design.csv > design.txt
% awk 'BEGIN { FS=","; OFS=" " } { $1=$1; print $0 }' contrasts.csv > contrasts.txt
% Text2Vest design.txt design.mat
% Text2Vest contrasts.txt design.con

% module load fsl
% https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/GLM/CreatingDesignMatricesByHand

%% Once the above is made, run the following in MATLAB

%palm -i data.csv -d design.mat -t design.con -o TermVsPreterm_TFCE -T -tfce2D -fdr -twotail -s week-40_hemi-left_space-dhcpSym_dens-32k_wm.surf.gii -logp 
%palm -i data.csv -d design.mat -t design.con -o TermVsPreterm_NO_TFCE -fdr -twotail -s week-40_hemi-left_space-dhcpSym_dens-32k_wm.surf.gii -logp
