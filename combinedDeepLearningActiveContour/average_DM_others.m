clc
clear all
close all
%%

% Ngo et al. 2014 
GC_a=[95.91,90.54];
GC_s=[5.28 14.4];

DM_a=[.88 .89];
DM_s=[.03, .03];

PD_a=[2.34,2.17];
PD_s=[.46,.46];

AGC=mean(GC_a);
sd_GC=mean(GC_s);
disp([num2str(AGC),'(',num2str(sd_GC),')'])

ADM=mean(DM_a);
sd_DM=mean(DM_s);
disp([num2str(ADM),'(',num2str(sd_DM),')'])


APD=mean(PD_a);
sd_PD=mean(PD_s);
disp([num2str(APD),'(',num2str(sd_PD),')'])


