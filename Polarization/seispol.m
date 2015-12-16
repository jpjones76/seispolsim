function varargout = seispol(X, varargin)
% P = seispol(X);
% [P, W] = seispol(X, av, Na, fca, z);
%
% Compute polarization attributes for X from Vidale (1986) + Jurkevics
% (1988). All attributes computed are returned in structure P. 
%
% INPUTS
% X     For K stations, 3K column vectors, arranged [Z_1 N_1 E_1 ...
%       Z_2 N_2 E_2 ... Z_K N_K E_K]. For Z-component data, down is +.
%
% OPTIONs   DEFAULT         MEANING
% av        1               Do time-averaging of az, el, in? (1 = yes)
% Na        5               Width of moving average filter in samples
% fca       0               Full circle az, inc? (1 = yes)
% z         1:1:180         Phase angles in degrees
%
% OUTPUTS
% P     Structure containing polarization information. 
%       STR     MEANING              RANGE | fca = 0     RANGE | fca = 1
%       az      Azimuth (a)         -90 <= az <= 90      -180 <= az <= 180
%       el      Ellipticity           0 <= el <= 1
%       in      Incidence angle (b) -90 <= in <= 90      -180 <= in <= 180
%       pl      Planarity             0 <= pl <= 1
%       rc      Rectilinearity        0 <= rc <= 1
%
% W     Weights for each measurement in P. 1x5 cell structure.
%
% NOTES ON OUTPUTS
% (a)   Possible abuse of terminology: azimuth is CLOCKWISE FROM NORTH, 
% thus +pi/6 = N30E, -pi/6 = N30W, etc. This has the same sense as strike;
% many polarization papers define azimuth as anticlockwise from East. 
%
% (b)   Incidence is FROM VERTICAL. This sense matches Jurkevics (1988),
% but is the complement of the "dip" parameter defined in Vidale (1986).  
%
% ======================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.5, 2015-12-16

%% Parameter initialization
av  = 1;                % Do time-averaging of az, el, in?
Na  = 5;                % Width of moving average filter
fca = 0;                % Full circle azimuth? 1 = yes
z   = 1:1:180;          % Vector of phase angles.

%% Parse varargin
if nargin > 1
    if ~isempty(varargin{1})
        av = varargin{1};
    end
    if nargin > 2
        if ~isempty(varargin{2})
            Na = varargin{2};
        end
        if nargin > 3
            if ~isempty(varargin{3})
                fca = varargin{3};
            end
            if nargin > 4 && ~isempty(varargin{4})
                z = varargin{4};
            end
        end
    end
end

[Nx,Nc] = size(X);
Nk = Nc/3;

%% Sanity checks
if rem(Nk,1); error('3C data required!'); end       % 3-component data only
if isreal(X); X = hilbert(X); end                   % Ensure X is analytic
Na = min([Na floor(Nx/3)]);                         % Filter order

%% Compute polarization attributes
P.az = zeros(Nx,1,Nk);
P.el = zeros(Nx,1,Nk);
P.in = zeros(Nx,1,Nk);
P.pl = zeros(Nx,1,Nk);
P.rc = zeros(Nx,1,Nk);
W1 = zeros(Nx,1,Nk);
W2 = zeros(Nx,1,Nk);
for k = 1:1:Nk
    x = X(:,3*k-2:3*k);
    C = [x(:,1).*conj(x(:,1)) x(:,1).*conj(x(:,2)) ...
        x(:,1).*conj(x(:,3)) x(:,2).*conj(x(:,1)) ...
        x(:,2).*conj(x(:,2)) x(:,2).*conj(x(:,3)) ...
        x(:,3).*conj(x(:,1)) x(:,3).*conj(x(:,2)) ...
        x(:,3).*conj(x(:,3))];
    
    % Jurkevics attributes
    C = circshift(filter((1/Na)*ones(1, Na), 1, C),[-floor(Na/2) 0]);
    CR = real(C);
    for m = 1:1:Nx
        [~,u] = eig([CR(m,1) CR(m,2) CR(m,3); ...
                     CR(m,4) CR(m,5) CR(m,6); ...
                     CR(m,7) CR(m,8) CR(m,9)],'nobalance');
        u = sort(diag(u),'descend');
        P.rc(m,k) = 1 - ( (u(2) + u(3)) / (2*u(1)) );
        P.pl(m,k) = 1 - ( 2*u(3) / (u(1)+u(2)) );
        W2(m,k) = C(1) + C(5) + C(9);
    end

    % Vidale attributes
    if av
        V1 = zeros(3,Nx);
        for m = 1:1:Nx-Na
            c = [C(m,1) C(m,2) C(m,3); ...
                 C(m,4) C(m,5) C(m,6); ...
                 C(m,7) C(m,8) C(m,9)];
            [V, d] = eig(c,'nobalance');
            [~, i1]  = sort(diag(d),'descend');
            v1 = V(:,i1(1));
            V1(:,m) = v1;
            W1(m,k) = sum(abs(diag(c)));
        end
        [P.az(:,k), P.el(:,k), P.in(:,k)] = seisvpol(V1, z, fca);
    else
        V = [x(:,1)./x(:,3) x(:,2)./X(:,3) ones(Nx,1)]';
        V = V./repmat(sqrt(sum(V.*conj(V))),3,1);
        [P.az(:,k), P.el(:,k), P.in(:,k)] = seisvpol(V, z, fca);
        W1(:,k) = sum(real(x).^2,2);
    end
end

%% Truncate 
Na = 1+ceil((Na-1)/2);
P.az = P.az(Na:Nx-Na,:);
P.el = P.el(Na:Nx-Na,:);
P.in = P.in(Na:Nx-Na,:);
P.pl = P.pl(Na:Nx-Na,:);
P.rc = P.rc(Na:Nx-Na,:);
W1 = W1(Na:Nx-Na,:);
W2 = W2(Na:Nx-Na,:);
P = orderfields(P);
W = cell(1,5);
W{1} = W1;
W{2} = W1;
W{3} = W1;
W{4} = W2;
W{5} = W2;

%% Parse varargout
varargout{1} = P;
if nargout > 1
    varargout{2} = W;
end

% ===============================================================
function [az, el, in] = seisvpol(V1, z, fca)
cis = (cosd(z) + 1i*sind(z))';
[xr, ir] = max(sqrt(real(cis*V1(1,:)).^2 + ...
                    real(cis*V1(2,:)).^2 + ...
                    real(cis*V1(3,:)).^2 ));
Xc = [V1(1,:).*cis(ir)' ; V1(2,:).*cis(ir)' ; V1(3,:).*cis(ir)' ];
Xr = real(Xc);
el = real(sqrt(1-xr.^2) ./ xr)';
if fca
    az = atan2d( Xr(3,:) , Xr(2,:) );                       % Strike
    in = atan2d( sqrt(Xr(2,:).^2 + Xr(3,:).^2) , Xr(1,:));  % Inc (90-dip)
else
    az = atand( Xr(3,:) ./ Xr(2,:) );                       % Strike
    in = atand( sqrt(Xr(2,:).^2 + Xr(3,:).^2) ./ Xr(1,:));  % Inc (90-dip)
end
