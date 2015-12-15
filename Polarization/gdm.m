function G = gdm(varargin)
% G = gdm(N,T,P);
% Form normalized ground distance matrix of inter-bin distances. Returns 
% a cell structure G with numel(T) NxN ground distances.
% 
% N    Number of bins.
% T    Threshold distance, in bins, below which ground distance 
%      is not computed. 
%      
%      SPECIAL CASES
%      T      Form of GDM
%      0      eye(Nk);
%      0<T<1  (1-b)*eye(NK) + b*ones(Nk);
%      1      ones(Nk);
% 
% P    Periodicity flag, size = size(T). (OPTIONAL) 
%      P = 1 if histograms correspond to a periodic function.
%      P = 0 otherwise. (Default)
%
% F    Bin distance function, cell struct of strings, size(T). (OPTIONAL)
%      'eye'  = Identity matrix, equivalent to chi-squared distance
%      'ones' = Matrix of 1s, ignore cross-bin distance completely
%      'const'= Matrix of 1s on diagonals and T(p) on off-diagonals
%      'lin'  = Linear falloff from 0 to T(p).
%      'gauss'= Normal distribution about 0. (Default)
%
% ====================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-12-14


N = varargin{1};
T = varargin{2};
P = zeros(size(T));
F = repmat({'gauss'},size(T));

if nargin > 2
    P = varargin{3};
    if nargin > 3
        F = varargin{4};
    end
end

Np = numel(T);
if numel(P) < Np
    P = [P repmat(P(end),[1 Np-numel(P)])];
end

G = cell(1,Np);
for p = 1:1:Np
    switch lower(F{p})
      case 'eye'
        G{p} = eye(N);
      case 'ones'
        G{p} = ones(N);
      case 'const'
        G{p} = T(p)*ones(N) + (1-T(p))*eye(N);
      case 'lin'
        if P(p)
            A = eye(N);
            G{p} = A;
            for a = 1:1:T-1
                G{p} = G{p} + ...
                       ( circshift(A,[0 -a]) + ...
                         circshift(A,[0 a]) ) .* ...
                       ((T(p)-a)/T(p));
            end
        else
            G{p} = zeros(N);
            for a = 1:1:N
                for b = max([1 a-T(p)+1]):1:min([N a+T(p)-1])
                    G{p}(a,b) = 1-(abs(a-b)/T(p)); 
                end
            end
        end
      
      % Default: Gaussian
      otherwise
        G{p} = zeros(N);
        if P(p)
            g = [normpdf(-T(p):1:T(p),0,T(p)/3) ...
                 zeros(1,N-2*T(p)-1)];
            g = g./max(g);
            g = circshift(g,[0 -T(p)]);
            for a = 1:1:N
                G{p}(a,:) = circshift(g,[0 a-1]);
            end
        else
            for a = 1:1:N
                b1 = max([1 a-T(p)+1]);
                b2 = min([N a+T(p)-1]);
                g = normpdf(b1:1:b2,a,floor(T(p)/3));
                G{p}(a,b1:b2) = g./max(g);
            end
        end
    end
end
