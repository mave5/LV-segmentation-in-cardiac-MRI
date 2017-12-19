%% Power Analysis plot
clc
clear all
close all
%%

for k=1:5
    L1=27*k;
    fn1=['matFiles/simResults/V_64_H100_TrainDataLength',num2str(L1),'.mat'];
    load (fn1,'dm_cv')

    DM(:,k)=dm_cv;
    
end
num_of_patients=[3,6,9,12,15]

x=[num_of_patients,18,21,24,27,30];
p1 =  -0.0002924;
p2 =0.01709;
p3 =0.692;
y = p1*x.^2 + p2*x + p3;
y = - 0.00029*x.^2 + 0.017*x + 0.69;

x=[1,2,3,4,5,6,7,8,9,10];

boxplot(DM,num_of_patients)
%hold on
%plot(x,y)

%axis([0 10 0.7 1])
set(gca,'XTick',0:3:30,'FontSize',16,...
   'FontName','Times New Roman');

mDM=median(DM);
hold on
plot(mDM,'+')
xlabel('Number of patients','FontSize',14)
ylabel('Dice Metric','FontSize',14)
title('Dice Metric by number of patients','FontSize',14)
box off

num_of_patients=[3 6 9 12 15];
p = polyfit(num_of_patients,mDM,2);
y=polyval(p,num_of_patients);
num_of_patients=[num_of_patients,18,21,24,27,30];
y=polyval(p,num_of_patients);
plot(y)
axis([0 10 0.7 1])
set(gca,'XTick',0:3:30,'FontSize',16,...
   'FontName','Times New Roman');
