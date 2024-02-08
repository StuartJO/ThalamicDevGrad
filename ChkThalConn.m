sub_ses = readtable('ALL_sub_ses.txt');

SUB = sub_ses.Var1;
SES = sub_ses.Var2;

for i = 1:length(SUB)
    
IN = load(['.\outputs\Weighted\',SUB{i},'_',num2str(SES(i)),'.mat'],'CortConnNorm','ThalConnNorm','ThalConn','CortConn');

CortConnNorm(i,:) = IN.CortConnNorm;
ThalConnNorm(i,:) = IN.ThalConnNorm;
ThalConn(i,:) = IN.ThalConn;
CortConn(i,:) = IN.CortConn;
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

Seeds100conn = mean(ThalConn)>100;
ThalConnMean = mean(ThalConn);
save('SeedsThr.mat','Seeds100conn','ThalConnMean')