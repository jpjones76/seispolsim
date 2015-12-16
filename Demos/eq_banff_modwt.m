% Test script of Dk, Dt using M = 3.1 SE of Banff
clear all; 
close all; 

disp('Read earthquake data');
evt = 'default';
set_hsim_params;

% Some default parameters are overridden to illustrate t-f behavior
X = read_test_evt(evt, fs, fc, 0);
setfigdefs;

% Take modwpt
[P, Wd, H, S, T] = polsim_modwt(X, J, wlt, fs, L, La);
h = pstf1(X, Wd, S, T, fs);
