function D = iwpt(W, pos, wlt, mo)
% function D = idwpt(w, [j n], wlt, mo);
%
% "Reverse" implementation of Wavelet Packet Transform. Creates the
% (MO)DW(P)T detail coefficients D{j,k} associated with a level j and band
% n.
%
% w	Input wavelet packet coefficient vector. w should be a row
%       vector, such as any one element of the cell structure W
%       produced by modwpt.m
%
% wlt	Wavelet to use (optional). Specify as e.g. 'D{2}' (Haar), 'D{4}',
% 	Daubechies wavelet of length 4, etc. Default is D{4}. Be sure that
% 	this wavelet exists in wletdb.mat before invoking it!
%
%
% m0    0 for inverse DWPT (default)
%       1 for inverse MODWPT
%
% ========================================================================
% Version: 1.2, last modified 2015-12-14
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com

J = pos(1);
N = pos(2);
if iscell(W)
    tmp = W{J,N};
    W = tmp;
    clear tmp;
end
[wf, ~, ~, C] = getwlt(wlt, J, mo);
c = C{J,N};
[Nx,Nk] = size(W);
if mo
    D = zeros(Nx,Nk);
else
    D = zeros(Nx*2^J, Nk);
end
for k = 1:1:Nk
    for j = J:-1:1
        if c(j) == 0
            f = wf{1};
        else
            f = wf{2};
        end
        if mo
            W1 = imodwptp(W(:,k), f, j);
        else
            W1 = idwptpyr(W(:,k), wf{1}, j);
        end
        W(:,k) = W1';
    end
    D(:,k) = W(:,k);
end
