function varargout = polhist(P, W, varargin)
% H = polhist(P);
% [H, D, T] = polhist(P, W, OPTIONS);
%
%        Compute histograms for polarization structure P, with weights W.
% 
% INPUTS
% P     Polarization structure with Np fields. Fields must be P.az, P.el,
%           P.in, P.rc, P.pl, ordered as shown.
%
% OUTPUTS
% H     Histograms for P, a 1 x Np cell structure.
% D     Histogram distances for P. See below for the structure.
% T     Times for histogram time slices (only valid with m = 't').
%
% OPTIONS
%       The following options can be specified. Always specify options as
% 'name', value pairs, e.g. polhist(P, 'm', 't', 'cs', [1 0 1 0 0])
%
% Name  Default         Type    Meaning
% m     '1'             str     Mode of histogram calculation:
%                                '1'    one histogram for each column in P
%                                't'    time slices of length L, advanced
%                                       La samples between calculations
% G     empty cell      cell    1xNp cell structure of ground dist.
%                                 matrices. Type 'help gdm' for details.
% L     100             int     Sample length for time slices
% La    20              int     Advance by La samples between time slices
% dp    0.01            num     Polarization histogram bin width,
%                                 normalized (0 < dp < 0.5)
% hd    0               bool    Compute distances? 1 = 'yes'
% W     ones            cell    Weights, size 1 x Np cell struct
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
% ========================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-12-14

%% Program defaults
m   = '1';
L   = 100;
La  = 20;
dp  = 0.01;
hd  = 0;
G   = cell(1,5);

%% Parse varargin
if numel(varargin) > 0
    j = 1;
    while j < numel(varargin)
        eval([varargin{j} '= varargin{j+1};']);
        j = j+2;
    end
end

% Ensure P, W are 3 dimensional
[P,W] = pexpand(P,W);

%% Call to subroutine
switch lower(m)
    case 't'
        [H, D, T]   = pht(P, W, G, hd, dp, L, La);
    otherwise
        [H, D]      = ph1(P, W, G, hd, dp);
        T           = [];
end

%% Parse varargout
varargout{1} = H;
if nargout > 1
    varargout{2} = D;
    if nargout > 2
        varargout{3} = T;
    end
end

% ================================================================
function [H, D, T] = pht(P, W, G, d, dp, L, La)
F           = fieldnames(P);
Np          = numel(F);             % # pol parameters
T           = [];                   % Time
x0          = 0:dp:1-dp;            % Left histogram border
x1          = [x0(2:end) 1];        % Right histogram border
D           = [];                   % Distances
[Nx, ~, Nk] = size(W{1});
Npk         = Np*Nk;
rx0         = repmat(x0,[1 1 Nk]);
rx1         = repmat(x1,[1 1 Nk]);

pov = (L-La)/L;
Nt = fix((Nx/L - pov) * (1/(1-pov)));
if d
    D = zeros(2*Npk, Nt);
end

% Compute H
for p = 1:1:Np
    pol = P.(F{p});
    j = 0;    
    for t = 0:La:Nx-L
        j = j+1;
        T(j) = t+1;
        H{p}(:,:,j) = squeeze( sum( ...
            repmat(W{p}(t+1:t+L,:,:),[1 1/dp 1]) .* ...
            bsxfun(@ge,pol(t+1:t+L,:,:),rx0) .* ...
            bsxfun(@lt,pol(t+1:t+L,:,:),rx1) ) );
        if d && t >= La
            h0 = squeeze( ...
                reshape(H{p}(:,:,[j-1 j]), 1/dp, 2*Nk, 1));
            d1 = hdist(h0,G{p});
            
            % Obtain D_t
            D(1 + (p-1)*Nk : p*Nk, j) = diag(d1(Nk+1:end,1:Nk));
            
            % All cross-station distances lag j by 1 unit; note change
            % of variable to time j-1 to compensate
            du = d1(1:Nk, 1:Nk);
            
            % Obtain average D_k between k, k1 ~= k
            D(Npk + 1 + (p-1)*Nk : Npk + p*Nk, j-1) = ...
                sum(du,2)./(Nk-1);
            
            % Obtain average D_k over all Nk*(Nk-1)/2 pairs
            D((2*Npk) + p, j-1) = ...
                mean(du(triu(du)~=0));
        end
    end
    if d
        D(1+(p-1)*Nk:p*Nk, 1) = D(1+(p-1)*Nk:p*Nk,2);
        
        % This copies 2nd-last distances to last distances.
        D(Npk + 1 + (p-1)*Nk : Npk + p*Nk, j) = ...
            D(Npk + 1 + (p-1)*Nk : Npk + p*Nk, j-1);
    end
end


% ================================================================
function [H, D] = ph1(P, W, G, d, dp)

F   = fieldnames(P);
Nk  = size(P.az, 3);
Np  = numel(F);
B   = 0:dp:1;
wh  = cell(Np);

H = cell(1,5);
if d; D = cell(1,5); end
for p = 1:1:Np
    pol = P.(F{p});
    wt = W{p};
    for k = 1:Nk
        Wn = wt(:,:,k);
        for j = 1:1:numel(B)-1
            i = pol(:,:,k) >= B(j) & ...
                pol(:,:,k) <  B(j+1);
            wh{p}(j,1) = sum(Wn(i));
            clear i
        end
        H{p}(:,k) = wh{p};
    end
   H{p}(isnan(H{p})) = 0;
    if d
        D{p} = hdist(H{p},G{p});
        D{p}(isnan(D{p})) = 0;
    end
end
