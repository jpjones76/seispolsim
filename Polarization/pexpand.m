function [P,W] = pexpand(P,W)
% P = pexpand(P)
% [P,W] = pexpand(P,W)
%
%   Expand all matrices in polarization structure P, weight structure W
% to three-dimensional arrays.
%
% ==============================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2015-12-14

F = fieldnames(P);
Np = numel(F);
for p = 1:1:Np
    if numel(size(P.(F{p}))) == 2
        pol = P.(F{p});
        [Nx, Nk] = size(pol);
        p2 = zeros(Nx, 1, Nk);
        for k = 1:1:Nk
            p2(:,1,k) = pol(:,k);
        end
        P.(F{p}) = p2;
    end
end

for p = 1:1:Np
    if numel(size(W{p})) == 2
        wt = W{p};
        [Nx, Nk] = size(wt);
        w2 = zeros(Nx, 1, Nk);
        for k = 1:1:Nk
            w2(:,1,k) = wt(:,k);
        end
        W{p} = w2;
    end
end
