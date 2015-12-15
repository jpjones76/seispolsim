function Xr = rotHoadley(varargin)
% Xr = rotHoadley(X);
% Xr = rotHoadley(X,S);
%
% Rotate data matrix X with seismograms from sensors S into empirical az, 
% inc. Reassign channel order to [Z N E] from [N' E' -Z].
%
% DEPENDENCIES: rotseis
%
% =======================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2014 09/29

X = varargin{1};
S = 1:1:12;
if nargin > 1
    S = varargin{2};
end

% Let Y = CSV data
Nc = size(X,2);
Nk = Nc/3;

% checked with seispol using 41 pt. averages
rotNew = [ ...
   -9.8857 ...
    4.4271 ...
  -50.3192 ...
  -10.7048 ...
  -28.7297 ...
   -0.7994 ...
  -16.7209 ...
  -24.8736 ...
   -6.0489 ...
   -1.5680 ...
  -29.1549 ...
   10.1344];

Y = zeros(size(X));
for k = 1:1:Nk
    Y(:,3*k-2) = -1*X(:,3*k);
    Y(:,3*k-1:3*k) = X(:,[3*k-2 3*k-1]);
end
Xr = rotseis(Y,rotNew(S));
