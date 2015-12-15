function [X,fs,st] = read_hoadley(varargin)
% [X,fs,st] = read_hoadley(fname);
%     Read Hoadley-specific CSV file fname and correct to true ground
% velocity. 
%
% [X,fs,st] = read_hoadley(fname,fs,fc);
%     Specify new sample rate (in Hz) and critical frequency (in Hz) for new
% instrument response. Defaults: fs = 100 Hz, fc = 3.0 Hz. fs must be an integer.
% 
% Outputs
% X = [z_1 n_1 e_1 z_2 ... n_k e_k];
% fs = sampling rate in Hz
% st = identifier string
%
% DEPENDENCIES: flatresp, rotHoadley

fs = 100;
f0 = 4000;
fc = 3;
fl = 1;
fh = 18;
fname = varargin{1};
if nargin > 1
    fs = varargin{2};
    if nargin > 2
        fc = varargin{3};
        if nargin > 3
            fl = varargin{4};
            if fl > 0
                fh = varargin{5};
            end
        end
    end
end

if ~strncmpi(fname(1),'/',1)
    fname = ['/data/Hoadley/Unsorted/' fname];
end
X = csvread(fname,8,1);

% Downsample, detrend
X = resample(detrend(X(:,1:36),'constant'),fs,f0);

% Rotate 
X = rotHoadley(X);

% Remove all columns that contain only zeros; note that searching for columns
% of zeros or using eps does not work due to weird data spec.
s = find(max(abs(X))<1.0e-10);

% Pull dead traces
if ~isempty(s)
    disp(['Removing columns ' num2str(s) ]);
    r = max(abs(X))>=1.0e-10;
    X = X(:,r);
end

% Flatten response
X = flatresp(X,fs,15,1,fc);
X = X*3.4169e-13; % Convert to m/s

% Filter
if fl
    disp(['Bandpassing data at ' num2str(fl) '-' num2str(fh) ' Hz']);
    [B,A] = butter(2,2*[fl/fs fh/fs]);
    X = filtfilt(B,A,X);
else
    disp('Not bandpassing; change by setting fl to a nonzero number');
end

[t1,r1] = strtok(fname,'/');
while ~isempty(r1)
    [t1,r1] = strtok(r1,'/');
end
st = regexprep(t1,'[ ,_.-]','');
st = strrep(st,'gdcma','');
