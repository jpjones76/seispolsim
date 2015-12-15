fs = 80;                            % Sampling frequency
wlt = 'LA{16}';                     % Wavelet filter to use
J = 5;                              % Max MODWT level
L =  fs;                            % Time window, samps
La = ceil(L/4);                     % Advance time, samps
fc = 3;                             % New fc for inst resp. Was fmin
fl = 1;                             % Low corner frequency
ctypes = {'hskt'};
target_dir_win = 'D:\Temp\';

% Initialize parameters for histogram scripts
dp = 0.01;
Nh = (1/dp);
thr = floor(Nh/6);
T = thr.*ones(1,5);
G = gdm(Nh,T,[1 0 1 0 0]);

