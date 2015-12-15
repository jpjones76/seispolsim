function tr = stalta(X, varargin)
% [s1,e1] = stalta(X,OPTS);
%
%    Very simple script for classical STA/LTA ratios of X. Initial LTA is
% set to mean(abs(X(1:L,:)). LTA is updated starting at t(k) > L whenever 
% mod(t(k-1),L) = 0.
%   If X has negative values or is complex, this script preprocesses X by
% taking the absolute value.
%
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-12-14

% Default options
Sw = 100;
Lw = 400;
dt = 1;
if nargin > 1
    Sw = varargin{1};
    if nargin > 2
        Lw = varargin{2};
        if nargin > 2
            dt = varargin{3};
        end
    end
end

if min(X) < 0
    X = abs(X);
end
[Nx,Nk] = size(X);
t = 0:dt:Nx-Sw;
T = numel(t);
L = mean(X(1:Lw,:),1);
A = Lw/dt;

% Output ratios in tr
tr = zeros(T,Nk);
for k = 1:1:T
    if t(k) > Lw
        if ~mod(k,A)
            L = mean(X(t(k-A)+1:t(k),:),1);
        end
    end
    tr(k,:) = mean(X(t(k)+1:t(k)+Sw,:),1) ./ L;
end
