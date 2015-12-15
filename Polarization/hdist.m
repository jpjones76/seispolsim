function D = hdist(H, G, varargin)
% D = hsim(H,G,OPTIONS);
%   Measure distance of histograms in matrix/structure H using
% metric t.
%
% REQUIRED INPUTS
% G     Ground distance matrix.
% H     Matrix or cell structure.
%       MATRIX: each column of H represents one histogram.
%       CELL:   each column of each cell of H represents
%               one histogram. All cells must contain same-sized
%               matrices (for now).
%
% OPTIONS
% Name  Default    Format        Notes
% t     'chi2'     str           Distance metric to use
%
% OUTPUTS
% D   Distances for H.
%     Format of distance matrix is
%     [D_p1k1 D_p1k2 ... D_p2k1 ... D_PK]
%     p is index to polarization quantity
%     k is index to station
%
%
% ============================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-12-13

t = 'qchd';
if numel(varargin) > 0
    t = varargin{1};
end

if iscell(H)
    [Nh,Nk] = size(H{1});
    M = numel(H);
    if M > 1
        h = H{1};
        for m = 2:M
            h(:,1+(m-1)*Nk:m*Nk) = H{m};
        end
        H = h;
        clear h;
    else
        H = H{1};
    end
else
    [Nh,Nk] = size(H);
    M = 1;
end

% Check for G
if ~iscell(G)
    tmp = G;
    clear G;
    G{1} = tmp;
    clear tmp;
end

D = [];
for m = 1:1:M
    d = zeros(Nk);
    h = H(:,1+(m-1)*Nk:m*Nk);
    
    switch lower(t)
        case 'chi2';
            % Chi-square distance (Snedcor & Cochran, 1967)
            h = h./repmat(sum(h),Nh,1);     % Normalize all h
            h(h<eps) = eps;                 % Prevent div/0
            for k1 = 1:Nh
                for k2 = k1+1:1:Nk
                    h1 = h(:,k1);
                    h2 = h(:,k2);
                    d(k1,k2) = 0.5*sum((h1-h2).^2 ./ (h1+h2));
                end
            end
            
        case 'qchd'
            % Quadratic Chi Histogram Distance (Pele & Werman, 2010)
            h = h./repmat(sum(h),Nh,1)/2  ;       % Normalize, max(D) ~ 1
            h(h<eps) = eps;                       % Zeros destablize
            
            % Compute QCHD
            for k1 = 1:Nh
                for k2 = k1+1:1:Nk
                    h1 = h(:,k1);
                    h2 = h(:,k2);
                    Z = (h1+h2)'*G{m};
                    Z(Z==0) = 1;
                    Z = Z.^0.5;                  % Can be 0 < pow < 1
                    hd = (h1-h2)'./Z;
                    d(k1,k2) = sqrt(max(hd*G{m}*hd',0));
                end
            end
            
        otherwise
            % Default to L2 norm
            for k1 = 1:Nh
                for k2 = k1+1:1:Nk
                    d(k1,k2) = norm(h(:,k1),h(:,k2));
                end
            end
    end
    D(:,1+(m-1)*Nk:m*Nk) = d + d';
end
