% Test script of Dk, Dt using M = 3.1 S of Banff
clear all; 
close all; 
evt = 'banff';
set_hsim_params;
X = read_test_evt(evt, fs, fc, fl);
setfigdefs;
figs = {'eq_banff_Sk.eps'; 'eq_banff_St.eps'};

% Labels for plots
Nk = size(X,2)/3;
lb = [1 2 4 5 6 7 8 9 10 12];
lbs = cell(1,Nk);
for k = 1:1:Nk
    lbs{k} = sprintf('%1s%02i','H',lb(k));
end

% Polarization similarity
[P,W] = seispol(X);
[H, D, T] = polhist(P,W,1,1);
S = adaptivesim(X, D, L, La);
h = psiplot(X, S, T, fs, L, [], lbs);