function varargout = modwt(X, wlt, J, varargin)
% W = modwt(X, wlt, J);
% [W, Lj, Sh] = modwt(X, wlt, J, wpt, pad);
%       Compute MODWT for X. 
% 
% If wpt = 1, compute full wavelet packet table.
% If pad = 1, pad X with zeros prior to computing modwt (prevents circular
%             boundary condition from affecting wavelet coeffs).
% 
% [W, Lj, Sh] = modwt(X, wlt, J);
%       Return additional outputs
% 
% Lj    Filter widths
% Sh    Circular shifts to make each W{j,n} zero-phase (1)
% 
% NOTES
% (1)   W{j,n} is only truly zero-phase with a Haar or Least Asymmetric
% wavelet filter ('symlets' in the wavelet toolbox). 
% 
% DEPENDENCIES: wscaledb.mat, modwptpy.c
%
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version 1.3, last modified 2015-12-14

wpt = 0;
pad = 0;
if nargin > 3
    wpt = varargin{1};
    if nargin > 4
        pad = varargin{2};
    end
end

% Ensure X is a column vector since we shell out to a .mex file.
[Nx, Nk] = size(X); 
if Nx < Nk;
    X = X';
    Nk = size(X,2);
end

% Get wavelet filters, etc
[wf, Lj, Sh, Seq] = getwlt(wlt, J, wpt);

% Pad X
if pad
    X = [X; zeros(Lj(J), Nk)];
end

% Initialize W
if wpt; N = 2^J; else N = 2; end;
W = cell(J, N);

% Construct wavelet packet table
for j = 1:J
    if wpt; N = 2^j; else N = 2; end;
    for n = 1:1:N
        if j > 1; X = W{j-1, ceil(n/2)}; end       
        
        % Filter with G or H?
        if Seq{j,n}(end) == 0
            filt = wf{1};
        else
            filt = wf{2};
        end
        
        % Call .mex
        for k = 1:1:Nk
            w = modwptpy(real(X(:,k)), filt, j);
            if ~isreal(X(:,k))
                wi = modwptpy(imag(X(:,k)), filt, j);
                w = w + 1i*wi;
            end
            W{j,n}(:,k) = w;
        end
    end
end

% parse varargout
if ~wpt;
    for j = 1:1:J-1
        W{j,1} = [];
    end
end
varargout{1} = W;
if nargout > 1
    varargout{2} = Lj;
    if nargout > 2
        varargout{3} = Sh;
    end
end
