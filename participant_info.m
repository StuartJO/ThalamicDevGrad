
load('UsedSubData.mat')

median(NonDup_birth_age)
min(NonDup_birth_age)
max(NonDup_birth_age)

%Nfemale = sum(strcmp(dHCP_meta.sex(in), 'female'));
Nfemale = sum(NonDup_sex==1);

median(NonTemplate_birth_age)
min(NonTemplate_birth_age)
max(NonTemplate_birth_age)

median(NonTemplate_scan_age)
min(NonTemplate_scan_age)
max(NonTemplate_scan_age)

sum(NonTemplate_sex==1)

sum(PretermTermScan_sex==1)
median(PretermTermScan_birth_age)
min(PretermTermScan_birth_age)
max(PretermTermScan_birth_age)

median(PretermTermScan_scan_age)
min(PretermTermScan_scan_age)
max(PretermTermScan_scan_age)