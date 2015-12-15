function varargout = xcalign(X, varargin)
% [Xa,Sh,C] = xcalign(X,L);
%
%    Cross-correlate the columns of X using an iterative least-squares
% approach, within allowed range of lags L. Return circularly shifted
% matrix that aligns the columns in time.
%
% INPUTS
% X    Input data with channels in column vectors.
% L    Lags (Optional). L can be a matrix or scalar:
%        Matrix: L(j,k) is the maximum allowed lag of X(:,j),X(:,k).
%        Scalar: L is the maximum allowed lag of any two columns.
%
% % [Xa,Sh,C] = xcalign(X,L,th,wm,nq);
%
% OPTIONAL INPUTS: Fine-tuning control
% th    Only use correlations whose values exceed th. DEFAULT: th = 0.0.
% wm    Weighting method.
%        'cd'   Weight based on correlation distance. Use if perfect
%        DEFAULT: w = max(xcov(X(:,j),X(:,k),lag)).
% nq    Number of autocovariance residuals for Ljung-Box Q test.
%       Only used if wm = 'lb'. DEFAULT: nq = ceil(log(length(X))).
%
% OUTPUTS
% Xa   Aligned data
% Sh   Matrix of shifts used to align each column in Xa.
%        Note: Sense of shifts is such that circshift(Xa(:,k),Sh(k))
%              aligns each column, not -Sh(k).
% C    Mean correlations of each (shifted) channel.
%
%
% REFERENCE
%     VanDecar and Crosson (1990), Determination of teleseismic
% relative phase arrival times using multi-channel cross-correlation
% and least squares, BSSA 80(1), 150-169.
%
% ====================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.2, 2012-11-26

thr = 0.0;
x_st = 1;
v = 1;
cdw = 0;

% Get X
X = detrend(X,'constant');

% Initialize counter to number of traces, P; ensure row vector.
[N,P] = size(X);
if P>N; X=X'; [N,P]=size(X); end

% Autoset L if unspecified
if nargin > 1
    L = varargin{1};
    if nargin > 2
        thr = varargin{2};
        if nargin > 3
            cdw = varargin{3};
        end
    end
else
    nq = ceil(log(N));
    L = N-nq-1;
end

% Truncate X if needed
if x_st > 1; X = X(x_st:end,:); N = N-x_st+1; end
% IMPORTANT NOTE: "Out of range" errors can result from x_st + L >
% length of original data. No way around this.

% Initialize counter to number of pairings
M1 = P*(P-1)/2;

% Set up M
M = zeros(M1,P);         % Start with zeros
M(M1+1,:) = 1;           % Constraint: dt sums to 0

% If a scalar was passed as L, expand L to a PxP matrix
if isscalar(L); L = L.*ones(P); end

% Set up dt
dt = zeros(M1,1);

% Set minimum weight and "seen" flag
eps1 = .01*thr;
seen = zeros(1,P);

% Fill M, dt, and weighting matrix
g = 0;
wt = zeros(1,M1);
for j = 1:P-1
    for k = j+1:P
        g = g+1;
        M(g,j) = 1;
        M(g,k) = -1;
        
        xc = xcov(X(:,j),X(:,k),'coeff');
        [xm,t] = max(xc(N-L(j,k)+1:N+L(j,k)-1));
        t = t-L(j,k);
        dt(g) = t;
        
        % Account for anticorrelation if cdw==1
        if cdw == 1
            xd = sqrt(xc.^2);
            xm = max(xd(N-L(j,k)+1:N+L(j,k)-1));
            wt(g) = xm;
        else
            wt(g) = xm;
        end
        
        if wt(g) >= thr
            seen(j) = 1;
            seen(k) = 1;
        end
    end
end

Xa = X;
if numel(find(seen==1)) >= (P-1)
    
    % Normalize wt
    wt = wt/max(wt);
    wt(wt<=eps1) = eps1;
    wt(g+1) = 1;
    wt = diag(wt);                  % Form diagonal matrix
    wt(isnan(wt)) = 1;
    dt(g+1) = 0;
    
    % Invert
    te = round((pinv(M'*wt*M))*M'*wt*dt);
    
    % Circularly shift X
    for j = 1:P
        Xa(:,j) = circshift(X(:,j),-te(j,1));
    end
else
    if v
        disp('Not aligning X. Inversion too poorly constrained');
    end
    te = zeros(1,P);
end

varargout{1} = Xa;
if nargout > 1
    varargout{2} = -te;
    if nargout > 2
        [C,~] = corrcoef(Xa);
        varargout{3} = sum(C - eye(size(C))) ./ (length(C)-1);
    end
end
