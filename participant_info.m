
[U,in] = unique(dHCP_meta.participant_id);
median(dHCP_meta.birth_age(in))
min(dHCP_meta.birth_age(in))
max(dHCP_meta.birth_age(in))

median(dHCP_meta.scan_age)
min(dHCP_meta.scan_age)
max(dHCP_meta.scan_age)

Nfemale = sum(strcmp(dHCP_meta.sex(in), 'female'));

full_birth_age = dHCP_meta.birth_age(I);
full_scan_age = dHCP_meta.scan_age(I);

[~,ORD] = sort(full_birth_age);

scatter(full_birth_age(ORD),1:length(I),'filled')
hold on
scatter(full_scan_age(ORD),1:length(I),'filled')