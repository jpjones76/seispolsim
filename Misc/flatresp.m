function varargout = flatresp(varargin)
% [Xt,r] = flatresp(X,fs,f0,h0,fc1,g0);
%
%      Translate X to a "flatter" instrument response curve, 
% based on fs, old rolloff frequency f0, and damping constant h.
%
% INPUTS
% X    Data
% fs   Sampling rate of X, Hz
% fc0  Original rolloff frequency, Hz
% h0   Original damping constant  (Optional)
% fc1  New rolloff frequency, Hz  (Optional)
% g0   Original gain              (Optional)
% 
% OUTPUTS
% Xt   Data translated to new critical frequency
% r    Instrument responses in CORAL-compatible format, given as 
%      [R_old R_new].
%      r(1)       Damping at rolloff
%      r(2)       Scalar gain in V/m/s
%      r(3)       # of complex poles (max. 28)
%      r(4:32)    complex poles
%      r(33)      # of complex zeros (max. 28)
%      r(34:62)   complex zeros
% 
% NOTES
% 1.) If f1 is not specified, f1 = 10^floor(log10(abs(f0/10))), 
%     start of next lowest decade. 
% 2.) Assumes a scaling constant 1/sqrt(2) at critical frequency.
% 3.) Instrument responses are structured for compatibility with 
%     CORAL by Ken Creager (kcc@ess.washington.edu). 
% 4.) If not specified, h = 1/sqrt(2).
% 
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2014-05-30

X = detrend(varargin{1},'constant');     % Remove mean on read
Fs = varargin{2};
f0 = varargin{3}; 
g = 1;
h = 1/sqrt(2);

% Use start of next lowest decade for new CF if f1 unspecified 
f1 = 10^floor(log10(abs(f0/10))); 

% Parse optional inputs
if nargin > 3
    h = varargin{4};
    if nargin > 4
        f1 = varargin{5};
        if nargin > 5
            g = varargin{6};
        end
    end
end

% Size and nextpow2 size of X
[Nx,Nc] = size(X);
N2 = 2^nextpow2(Nx);  

% Initialize response vectors
R0 = zeros(62,1);
R1 = R0;

% Old instrument response
R0(1) = h;
R0(2) = g;
R0(3) = 2;
R0(4) = (-h + sqrt(h^2-1))*2*pi*f0;
R0(5) = (-h - sqrt(h^2-1))*2*pi*f0;
R0(33) = 2;
R0(34:35) = 0;

% New instrument response
h1 = 1/sqrt(2);
R1(1) = h1;
R1(2) = g;
R1(3) = 2;
R1(4) = (-h1 + sqrt(h1^2-1))*2*pi*f1;
R1(5) = (-h1 - sqrt(h1^2-1))*2*pi*f1;
R1(33) = 2;
R1(34:35) = 0;

% Take FFT at frequencies in f
f = [0:N2/2   -N2/2+1:-1]'/(N2/Fs);
Xf = fft(X,N2);

% New response at f
Rn = freqs( ...
    poly(R1(34:R1(33)+33)), ...
    poly(R1(4:R1(3)+3)), ...
    2*pi*f)*R1(1);

% Old response at f
Ro = freqs( ...
    poly(R0(34:R0(33)+33)), ...
    poly(R0(4:R0(3)+3)), ...
    2*pi*f)*R0(1);

% Translate response
wl = 1.e-6;                   % Waterlevel to prevent div/0 
temp1 = Ro.*conj(Ro);
gamma = max(temp1)*wl;
temp = Rn.*conj(Ro)./(temp1+gamma);
Xf = real(ifft(Xf.*repmat(temp,[1 Nc])));
Xt = Xf(1:Nx,:);

% Correct gain and renormalize
Xt = Xt/R0(2);

% All done.
if nargout
    varargout{1} = Xt;
    if nargout > 1
        varargout{2} = [R0 R1];
    end
end
