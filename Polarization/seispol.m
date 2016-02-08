function varargout = seispol(X, varargin)
% P = seispol(X);
% [P, W] = seispol(X, av, Na, z);
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
% Na        13              Width of moving average filter in samples
% z         1:1:180         Phase angles in degrees
%
% OUTPUTS
% P     Structure containing polarization information. 
%       STR     MEANING              RANGE
%       az      Azimuth (a)         -90 <= az <= 90
%       el      Ellipticity           0 <= el <= 1
%       in      Incidence angle (b)   0 <= in <= 90
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
% Version: 1.6, 2016-1-15

%% Parameter initialization
av  = 1;                % Do time-averaging of az, el, in?
Na  = 21;               % Width of moving average filter
z   = 1:1:180;          % Vector of phase angles.

%% Parse varargin
if numel(varargin) > 0
    if ~isempty(varargin{1})
        av = varargin{1};
    end
    if nargin > 2
        if ~isempty(varargin{2})
            Na = varargin{2};
        end
        if nargin > 3 && ~isempty(varargin{4})
            z = varargin{4};
        end
    end
end

cis = (cosd(z) + 1i*sind(z))';
[Nx,Nc] = size(X);
Nk = Nc/3;

%% Sanity checks
if rem(Nk,1); error('3C data required!'); end       % 3-component data only
if isreal(X); X = hilbert(X); end                   % Ensure X is analytic
Na = min([Na floor(Nx/3)]);                         % Filter order

%% Compute polarization attributes
P.az = zeros(Nx,Nk);
P.el = zeros(Nx,Nk);
P.in = zeros(Nx,Nk);
P.pl = zeros(Nx,Nk);
P.rc = zeros(Nx,Nk);
W1 = zeros(Nx,Nk);
W2 = zeros(Nx,Nk);
V = complex(zeros(3,Nx));
U = zeros(Nc,Nx);

for k = 1:1:Nk
    x = X(:,3*k-2:3*k);
    C = [repmat(x(:,1),1,3).*conj(x(:,1:3)) ...
         repmat(x(:,2),1,3).*conj(x(:,1:3)) ...
         repmat(x(:,3),1,3).*conj(x(:,1:3))];
    
    % Jurkevics attributes
    C = circshift(filter((1/Na)*ones(1, Na), 1, C),[-floor(Na/2) 0]);
    for m = 1:1:Nx
        c = [C(m,1:3); C(m,4:6); C(m,7:9)];
        [v, d] = eig(c,'nobalance');
        [d, i]  = sort(diag(d),'descend');
        if av
            V(:,m) = v(:,i(1));
        end
        U(3*k-2:3*k,m) = d;
        W2(m,k) = sum(diag(c));
    end

    % Vidale attributes
    if ~av
        V = [x(:,1)./x(:,3) x(:,2)./x(:,3) ones(Nx,1)]';
        W1(:,k) = sum((x.*conj(x)),2);
        V = V./repmat(sqrt(sum(V.*conj(V),1)),3,1);
    end
    [P.az(:,k), P.el(:,k), P.in(:,k)] = seisvpol(V, cis);
end
U = U';
U = U(1:Nx-Na,:);
P.rc = (1 - ( (U(:,2:3:end) + U(:,3:3:end)) ./ (2*U(:,1:3:end)) ));
P.pl = (1 - ( 2.*U(:,3:3:end) ./ (U(:,1:3:end)+U(:,2:3:end)) ));

%% Truncate 
P.az = P.az(1:Nx-Na,:);
P.el = P.el(1:Nx-Na,:);
P.in = P.in(1:Nx-Na,:);
W1 = W1(1:Nx-Na,:);
W2 = W2(1:Nx-Na,:);
P = orderfields(P);
W = cell(1,5);
if av
    W{1} = W2;
    W{2} = W2;
    W{3} = W2;
else
    W{1} = W1;
    W{2} = W1;
    W{3} = W1;
end
W{4} = W2;
W{5} = W2;

%% Parse varargout
varargout{1} = P;
if nargout > 1
    varargout{2} = W;
end

% ===============================================================
function [az, el, in] = seisvpol(V, cis)
[xr, ir] = max(real(cis*V(1,:)).^2 + real(cis*V(2,:)).^2 + ...
               real(cis*V(3,:)).^2, [], 1);
el = real(sqrt((1-xr)./xr))';
Xr = real(V .* repmat(cis(ir)',3,1))';
az = atand( Xr(:,3) ./ Xr(:,2) );
in = abs(atand(sqrt(Xr(:,2).^2 + Xr(:,3).^2) ./ Xr(:,1)));
