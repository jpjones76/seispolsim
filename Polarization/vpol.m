function varargout = vpol(X, varargin)
% [az, el, in] = vpol(X, fca, z);
%   Standalone to compute Vidale (1986) instantaneous polarization
% attributes (azimuth, ellipticity, and incidence). Returns column vectors,
% one column per station.
%
% [az, el, in, W] = vpol(X, z, fca);
%   Also return the summed energies of X as weights W
%
% REQUIRED INPUT
% X     3c seismic data arranged [Z_1 N_1 E_1 Z_2 N_2 E_2 ... Z_K N_K E_K]
%
% OPTIONAL INPUTS
% fca   Compute full-circle angular attributes? (1 = yes, default = 0)
% z     Row vector of complex phase angles in degrees (default = 1:1:180)
%
% OUTPUTS           RANGE
%                   fca = 0     fca = 1
% el   Ellipticity  [0, 1]      [0, 1]
% az   Azimuth      [-90, 90]   [-180, 180]
% in   Incidence    [-90, 90]   [-180, 180]
% 
% Standalone
%
% ======================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, 2015-12-10

z = 1:1:180;
fca = 0;
if numel(varargin) > 0
    z = varargin{1};
    if numel(varargin) > 1
        fca = varargin{2};
    end
end

% Form normalized X
[Nx, Nc] = size(X);
if rem(Nc,3)
    error([mfilename ' requires 3-component data!']);
else
    Nk = Nc/3;
end

% Analytic signal
if isreal(X)
    X = hilbert(X);
end

W = zeros(Nx, 1, Nk);
el = zeros(Nx, Nk);
az = zeros(Nx, Nk);
in = zeros(Nx, Nk);
cis = (cosd(z) + 1i*sind(z))';

for k = 1:1:Nk
    V = [X(:,3*k-2)./X(:,3*k) X(:,3*k-1)./X(:,3*k) ones(Nx,1)]';
    V = V./repmat(sqrt(sum(V.*conj(V))),3,1);
    
    % Eqn. (4) of Vidale (1986), including rotation through phase angle
    [xr, ir] = max( sqrt( real(cis*V(1,:)).^2 + real(cis*V(2,:)).^2 + ...
                          real(cis*V(3,:)).^2 ) );
    Xc = [V(1,:).*cis(ir)' ; ...
          V(2,:).*cis(ir)' ; ...
          V(3,:).*cis(ir)' ];
    Xr = real(Xc);

    % Ellipticity
    el(:,k) = real(sqrt(1-xr.^2)./xr)';
    
    % Azimuth and Incidence
    if fca
        az(:,k) = atan2d( Xr(3,:) , Xr(2,:) );
        in(:,k) = atan2d( Xr(1,:) , sqrt(Xr(2,:).^2 + Xr(3,:).^2) );
    else
        az(:,k) = atand( Xr(3,:) ./ Xr(2,:) );
        in(:,k) = atand( Xr(1,:) ./ sqrt(Xr(2,:).^2 + Xr(3,:).^2) );
    end
end

varargout{1} = az;
varargout{2} = el;
varargout{3} = in;
if nargout > 3
    for k = 1:1:Nk
        W(:,1,k) = sum(real(X(:,3*k-2:3*k).^2),2);
    end
    varargout{4} = W;
end
