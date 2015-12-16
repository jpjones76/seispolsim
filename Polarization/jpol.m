function [az, in, pl, rc, W] = jpol(X, N, varargin)
% [az, in, pl, rc] = jpol(X, N);
%
%   Compute Jurkevics' time-averaged polarization attributes. Return
% planarity and rectilinearity as column vectors, one column per station.
% An N-point moving average filter is used to stabilize the calculations.
%
% [az, in, pl, rc] = jpol(X, N, fca);
%   Set fca to 1 to resolve 180 degree ambiguity in az, in using atan2d
% 
% [az, in, pl, rc, W] = vpol(X, z, fca);
%   Also return summed energies of X as weights W
%
% REQUIRED INPUT
% X     3c seismic data arranged [Z_1 N_1 E_1 Z_2 N_2 E_2 ... Z_K N_K E_K]
%
% OUTPUTS               RANGE
% az   Azimuth          [-180, 180] (fca = 1)
%                       [-90, 90]   (fca = 0)
% in   Incidence        [0, 90]
% rc   Rectilinearity   [0, 1]
% pl   Planarity        [0, 1]
%
% % Standalone
%
% ======================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-12-10

[Nx, Nc] = size(X);
Nk = Nc/3;
os  = ceil((N-1)/2);
az = zeros(Nx, Nk);
in = zeros(Nx, Nk);
rc = zeros(Nx, Nk);
pl = zeros(Nx, Nk);
fca = 0;
if numel(varargin) > 0
    fca = varargin{1};
end
if ~isreal(X)
    X = real(X);
end

W = zeros(Nx, Nk);
for k = 1:1:Nk
    x = X(:,3*k-2:3*k);
    C = [x(:,1).^2      x(:,1).*x(:,2) x(:,1).*x(:,3) ...
         x(:,2).*x(:,1) x(:,2).^2      x(:,2).*x(:,3) ...
         x(:,3).*x(:,1) x(:,3).*x(:,2) x(:,3).^2];
    C = filtfilt(repmat(1/N,1,N),1,C);
    for m = 1+os:1:Nx-os
        [V,u] = eig([C(m,1) C(m,2) C(m,3); ...
                     C(m,4) C(m,5) C(m,6); ...
                     C(m,7) C(m,8) C(m,9)],'nobalance');
        [u,i1] = sort(diag(u),'descend');
        V = V(:,i1);
        v1 = V(:,1);
        
        % Polarizations are computed from first eigenvector
        if fca
            az(m,k) = atan2d( v1(2), v1(3) ); % No need to use signum
        else
            az(m,k) = atand( v1(2) / v1(3) );
        end
        in(m,k) = acosd( abs(v1(1)) );
        rc(m,k) = 1 - ( (u(2) + u(3)) / (2*u(1)) );
        pl(m,k) = 1 - ( 2*u(3) / (u(1)+u(2)) );
        
        % Weight
        if m > os && m <= Nx-os
            W(m-os,1,k) = C(1) + C(5) + C(9);
        end
    end
end
