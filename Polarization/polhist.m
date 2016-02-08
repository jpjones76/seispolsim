function varargout = polhist(P, W, varargin)
% H = polhist(P,W);
% [H, D, T] = polhist(P, W, ts, hd, L, La, G);
% [H, D, T] = polhist(P, W, ts, hd, L, La, Nb);
%
%        Compute histograms for polarization structure P, with weights W.
%
% INPUTS
% P     Polarization structure with Np fields. Fields must be P.az, P.el,
%           P.in, P.rc, P.pl, ordered as shown.
%
% OUTPUTS
% H     Energy-weighted histograms of the polarization of X, divided into
%       bins of normalized width dp. Order matches alphabetical order of
%       field in P: H{1} = az, H{2} = el, H{3} = in, H{4} = pl, H{5} = rc.
%
% D     Histogram distances. For Nk stations (3Nk channels), Nt time
%       slices, Np attributes, the distance matrix D is size Np(2Nk+1) x Nt.
%   ROW #       Meaning
%   1:NkNp          Histogram distances for attribute p, between time
%               slices i, i-1 (i = 2:Nt). Row order:
%                           p=1, k=1
%                           p=1, k=2
%                               ...
%                           p=Np, k=Nk-1
%                           p=Np, k=Nk
%                   
%   NkNp+1:         Mean histogram distance between histogram at time t,
%    2NkNp     station k, and histograms at time t, stations j != k, for
%               each attribute p
%                           p=1, k=1        (vs. p=1, k=2:Nk)
%                           p=1, k=2        (vs. p=1, k=1, 3:Nk)
%                               ...
%                           p=Np, k=Nk-1    (vs. p=Np, k=1:Nk-2, k=Nk)
%                           p=Np, k=Nk      (vs. p=Np, k=1:Nk-1)
%  2NkNp+1:
%   Np(2Nk+1)       Mean of all distances between histogram pairs of p at
%               time t
%
% T     Starting sample number of each histogram. Empty unless ht = .%
%
% OPTIONS
%       The following options can be specified. 
%
% Name  Default         Type    Meaning
% ts    0               bool     Mode of histogram calculation:
%                                0      One histogram for each P{p}(:,k)
%                                1      Histograms of time slices
% sim   0               bool    Compute histogram similarity? (1 = yes)
% G     empty cell      cell    1xNp cell structure of ground dist.
%                                 matrices. Type 'help gdm' for details.
% L     100             int     Sample length for time slices
% La    20              int     Advance by La samples between time slices
% Nb    100             int     Bins per histogram
%
% STRUCTURE OF D
%   For K stations (3K channels), the distance matrix D is Np(2K+1) x Nt,
% where Np is # of polarization attributes computed (either 3 or 5) and Nt
% is # of time slices computed.
%
%   ROW #       Meaning / Ordering
%   1:KNp       Histogram distance between time slice t and slice t-1 at
%               each station k for each attribute p
%                   ORDERING:   r[1-Np]: p=1, k=1:K; ...
%                               r[(Np-1)*Nk+1:Np*Nk] =  p=Np, k=1:K
%   KNp+1:2KNp  Mean histogram distance between histogram at time t,
%               station k, and histograms at time t, stations j != k, for
%               each attribute p
%                   ORDERING:   r[1-Np]: p=1, k=1:K; ...
%                               r[(Np-1)*Nk+1:Np*Nk] =  p=Np, k=1:K
%   2KNp+1:     Average Dk over all K station pairs for each attribute p
%     2KNp+Np       ORDERING: p = 1:Np
% 
% Dependencies: pexpand, gdm, hdist
% 
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.2, 2015-12-16

%% Program defaults
ts  = 0;
hd  = 0;
L   = 100;
La  = 20;
Nb  = 100;

%% Parse varargin
if numel(varargin) > 0
    if ~isempty(varargin{1})
        ts = varargin{1};
    end
    if numel(varargin) > 1
        if ~isempty(varargin{2})
            hd = varargin{2};
        end
        if numel(varargin) > 2
            if ~isempty(varargin{3})
                L = varargin{3};
            end
            if numel(varargin) > 3
                if ~isempty(varargin{4})
                    La = varargin{4};
                end
                if numel(varargin) > 4
                    if iscell(varargin{5})
                        G = varargin{5};
                        Nb = size(G{1},1);
                    elseif ~isempty(varargin{5})
                        Nb = varargin{5};
                    end
                end
            end
        end
    end
end

% Ensure P, W are 3 dimensional
[P,W] = pexpand(P,W);

% Dependent constants, arrays
[Nx, ~, Nk] = size(W{1});
F   = fieldnames(P);
Np  = numel(F);
Npk = Np*Nk;
T   = 1+(0:La:Nx-L);
Nt  = numel(T);
dp  = 1/Nb;
% pov = (L-La)/L;
% Nt = fix((Nx/L - pov) * (1/(1-pov)));

if hd && ~exist('G','var')
    G = gengdm(Nb);
end

% Histogram boundaries
x0          = 0:dp:1-dp;
x1          = [x0(2:end) 1];
rx0         = repmat(x0,[1 1 Nk]);
rx1         = repmat(x1,[1 1 Nk]);

% Initialization
H = cell(1,Np);
if hd && ts
    D = zeros(2*Npk+Np, Nt);
elseif hd
    D = zeros(2*Npk,Nk);
else
    D = [];
end

% Map az, inc to [0, 1]
P.az = P.az/180 + 0.5;
P.in = P.in/90;

% Compute H, D
for p = 1:1:Np
    pol = P.(F{p});
    pol(pol==1) = 1-eps;
    if ts
        j = 0;
        for t = 0:La:Nx-L
            j = j+1;
            H{p}(:,:,j) = squeeze( sum( ...
                repmat(W{p}(t+1:t+L,:,:),[1 1/dp 1]) .* ...
                bsxfun(@ge,pol(t+1:t+L,:,:),rx0) .* ...
                bsxfun(@lt,pol(t+1:t+L,:,:),rx1) ) );
            if hd && t >= La
                h0 = squeeze( ...
                    reshape(H{p}(:,:,[j-1 j]), 1/dp, 2*Nk, 1));
                d1 = hdist(h0, G{p});
                
                % Obtain D_t
                D(1 + (p-1)*Nk : p*Nk, j) = diag(d1(Nk+1:end,1:Nk));
                
                % Obtain average D_k for k1 ~= k
                du = d1(1:Nk, 1:Nk);
                D(1+Npk+(p-1)*Nk : Npk+p*Nk, j-1) = sum(du,2)./(Nk-1);
                
                % Obtain average D_k over all Nk*(Nk-1)/2 pairs
                D((2*Npk) + p, j-1) = mean(du(triu(du)~=0));
            end
        end
    else
        H{p} = squeeze(sum(repmat(W{p},[1 1/dp 1]).* ...
                       bsxfun(@ge,pol,rx0) .* ...
                       bsxfun(@lt,pol,rx1),1));
        H{p}(isnan(H{p})) = 0;
        if hd
            D(1+Npk+(p-1)*Nk : Npk+p*Nk, :) = hdist(H{p}, G{p});
        end
    end
end

% Check for single-station case with ts = 1 (squeeze is inconsistent)
if Nk == 1 & ts
    for p = 1:1:5
        H{p} = permute(H{p},[2 1 3]);
    end
end

% Remap az, in
P.az = 180*(P.az-0.5);
P.in = P.in*90;

%% Parse varargout
varargout{1} = H;
if nargout > 1
    varargout{2} = D;
    if nargout > 2
        varargout{3} = T;
    end
end

% ========================================================================
function D = hdist(H, G)
% D = hdist(H,G,OPTIONS);
%   Measure qchd of histograms in matrix/structure H.
%
% INPUTS
% G     Ground distance matrix.
% H     Matrix of histograms in column vectors.
%
% OUTPUTS
% D   Distances between all pairs in H.
%
% ============================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.2, 2016-01-28

[Nb, Nh] = size(H);
D = zeros(Nh);
H = 0.5*H./repmat(sum(H,1),Nb,1);
H(H<eps) = eps;
i0 = 0;
for j = 1:1:Nh-1
    i1 = i0+Nh-j;
    D(j+1:Nh,j) = qchd(repmat(H(:,j),1,Nh-j), H(:,j+1:Nh), G);
    i0 = i1;
end
D = D + D';