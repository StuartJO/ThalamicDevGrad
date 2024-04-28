load('UsedSubData.mat')

for i = 1:length(SUB)

data = dlmread(['D:/TC_connectivity/',SUB{i},'_',num2str(SES(i)),'_thal_conn_verts_wei.txt']);
data_nonmed = data(:,medwallmask);
data_nonmed(isnan(data_nonmed)) = 0;
Norm = BF_NormalizeMatrix(data_nonmed,'scaledSigmoid');
Norm(isnan(Norm)) = 0;

CortConnNorm(i,:) = sum(Norm,1);
ThalConnNorm(i,:) = sum(Norm,2);

ThalConn(i,:) = sum(data_nonmed,2);
CortConn(i,:) = sum(data_nonmed,1);

end

vox_coords = dlmread('thal_seed_1.75mm_vox_coords.txt');

ThalConnNormMean = mean(ThalConnNorm);

for i = 1:3
    subplot(1,3,i)
    scatter(ThalConnNormMean,vox_coords(:,i))
end

for i = 1:3
disp([num2str(mean(AllExpl(Test,i))),'+-',num2str(std(AllExpl(Test,i)))])
end

Seeds100conn = mean(ThalConn)>=100;
ThalConnMean = mean(ThalConn);
save('SeedsThr.mat','Seeds100conn','ThalConnMean')