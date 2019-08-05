% plot correlation and Bland Altman plots

clc
close all
clear all
%%

X=xlsread('res1.xlsx');
%X=xlsread('res2.xlsx');

% manual
EF_man=X(1:15,13);
EDV_man=X(1:15,23);
ESV_man=X(1:15,27);

% automatic
EF_auto=X(1:15,9);
EDV_auto=X(1:15,21);
ESV_auto=X(1:15,25);

% Bland-Altman
BlandAltman(EF_man,EF_auto,{'Manual EF(%)','Auto EF(%)'},'Ejection Fraction');
BlandAltman(EDV_man,EDV_auto,{'Manual EDV(cm^3)','Auto EDV(cm^3)'},'End-diastolic Volume');
BlandAltman(ESV_man,ESV_auto,{'Manual ESV(cm^3)','Auto ESV(cm^3)'},'End-systolic Volume');



