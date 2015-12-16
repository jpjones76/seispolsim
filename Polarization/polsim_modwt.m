function varargout = polsim_modwt(X, J, wlt, fs, varargin)
% D = polsim_modwt(X,J,wlt,dp);
%    MODWT for polarization similarity. Returns distance matrix D.
% 
% [D,P,H,W,Wd] = polsim_modwt(X,J,wlt,dp);
%    Full list of possible outputs. 
%
% INPUTS
% X    Data arranged in column vectors as [Z1 N1 E1 ... ZK NK EK].
% J    Maximum wavelet transform level J.
% wlt  Wavelet to use, syntax 'LA{8}'.
% dp   Polarization histogram spacing, e.g. 0.02.
%
% OUTPUTS
% D    Histogram distances
% P    Polarization structure
% H    Histogram structure, order matches that of P
% W    MODWT structure with wavelet coefficients in W{j,2} (j=1:J) and
%      scaling coefficients in W{J,1}
% Wd   Wavelet detail coefficients with same structural arrangement as W.
%
% DEPENDENCIES: getwlt, modwt, modwtpy, iwpt, seispol, wscaledb, polhist, 
% adaptivesim, stalta, imodwptp
%
% ==============================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-12-14

L = fs;
La = ceil(L/2);
if numel(varargin) > 0
    L = varargin{1};
    if numel(varargin) > 1
        La = varargin{2};
    end
end

Nc = size(X,2);
[~, Lj] = getwlt(wlt, 1, J);
X = [X; zeros(max(Lj),Nc)];

% Take modwt with X padded
W = modwt(X, wlt, J, 0, 1);

% Initialize cell structures
Wd = cell(J,2);
P = cell(J,2); 
H = cell(J,2);
D = cell(J,2);

% Do inverse modwt 
for j = 1:1:J; 
    if j == J
        nmin = 1;
    else
        nmin = 2;
    end
    for n = 2:-1:nmin
        Wd{j,n} = iwpt(W{j,n}, [j n], wlt, 1);
        [P{j,n}, wgt] = seispol(Wd{j,n});
        [~, D, T] = polhist(P{j,n}, wgt, 1, 1);
        S{j,n} = adaptivesim(Wd{j,n}, D, L, La);
    end; 
end


%% Parse varargout
varargout{1} = P;
if nargout > 1
    varargout{2} = Wd;
        varargout{3} = H;
        if nargout > 3
            varargout{4} = S;
            if nargout > 4
                varargout{5} = T;
            end
        end
    end
end