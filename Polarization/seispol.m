function varargout = seispol(X, varargin)
% P = seispol(X,OPTS);
% [P, W, H, D, T] = seispol(X,OPTS);
%
% Compute polarization for X by the method of Vidale(1986) +
% Jurkevics(1988). Attributes computed are returned in P.
%
% INPUTS
% X      [Z_1 N_1 E_1 Z_2 N_2 E_2 ... Z_K N_K E_K];
%
%       For K stations, 3K column vectors, arranged as above. For
% Z-component data, down is +.
%
%
% OUTPUTS
% P     Structure containing polarization information:
%       VAR    MEANING              RANGE | fca = 0     RANGE | fca = 1
%       P.az   Azimuth (a)         -90 <= az <= 90      -180 <= az <= 180
%       P.el   Ellipticity           0 <= el <= 1
%       P.in   Incidence angle (b) -90 <= in <= 90      -180 <= in <= 180
%       P.pl   Planarity             0 <= pl <= 1
%       P.rc   Rectilinearity        0 <= rc <= 1
%
% W     Weights for each measurement in P. 1x5 cell structure.
%
% H     Energy-weighted histograms of the polarization of X, divided into
%       bins of normalized width dp. Order matches alphabetical order of
%       field in P: H{1} = az, H{2} = el, H{3} = in, H{4} = pl, H{5} = rc.
%
% D     Distances between histograms. Empty unless hd = 1.
%
% T     Middle sample number of each histogram. Empty unless ht = .
%
% NOTES ON OUTPUTS
% (a)   I intentionally commit an abuse of terminology here: azimuth is
% measured CLOCKWISE FROM NORTH, thus +pi/6 = N30E, -pi/6 = N30W, etc. This
% is a strike direction; most polarization papers treat azimuth as
% anticlockwise from East. 
%
% (b)   Incidence is measured FROM VERTICAL, thus 0 = vertical. This sense
% matches Jurkevics (1988), but is the complement of the "dip" parameter
% defined in Vidale (1986).
%
%
% OPTS  (DEFAULT)         Meaning
% z   = 1:1:180           % Phase angles in degrees
% L   = 100;              % # of samples per hist segment
% La  = 20;               % # of samples to advance per segment
% Na  = 5;                % Width of moving average filter
% av  = 1;                % Do time-averaging of az, el, in?
% dp  = 0.01;             % Bin spacing for histograms
% ew  = 1;                % Energy-weight histograms? (Recommended)
% fca = 0;                % Full circle az, inc? 1 = yes
% h1  = 0;                % Do one-shot histograms? 1 = yes
% hd  = 0;                % Do histogram distances?
% ht  = 0;                % Do histogram time slices? 1 = yes
% pv  = 16;               % Variance of Gaussian used in GDM
% v   = 0;                % Verbosity
%
% TIPS FOR EFFECTIVE USE
% 1) Setting 'h1',1 and 'ht',1 computes "time slice" histograms only; this
% is exactly the same behavior as 'h1',1,'ht',0. (Only one output field for
% histograms)
%
% 2) If using repeat calls to seispol to construct "time slice" histograms
% for successive long windows of length Nx1 < Nx, use an overlapping data
% window. If La = L, advance time at each successive call to seispol by
%    t(end) - Na/2 + Nsamp/2
% rather than by Nx1 samples.
%
% CHANGELOG
% 2015 03-29    Incidence angle is now 0 <= in <= pi/2 if fca = 0
% 2015 04-01    Bugfix for h1 = 1
% 2015 04-15    Added info. to this help file.
% 2015 04-19    Moved histograms to external functions; code cleanup/speed
% 2015 04-20    Changed order of varargout
%
% ======================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.4, 2015-12-15

%% Parameter initialization
z   = 1:1:180;          % Vector of phase angles.
L   = 100;              % # of samples per hist segment
La  = 20;               % Sample advance
Na  = 5;                % Width of moving average filter
av  = 1;                % Do time-averaging of az, el, in?
dp  = 0.01;             % Bin spacing for histograms
fca = 0;                % Full circle azimuth? 1 = yes
h1  = 0;                % Do one-shot histograms?
ht  = 0;                % Do histogram time slices?
hd  = 0;                % Do histogram distances?
pv  = floor(1/(6*dp));  % Variance of Gaussian used in GDM

%% Parse varargin
if nargin > 1
    j = 1;
    while j < numel(varargin)
        eval([varargin{j} '= varargin{j+1};']);
        j = j+2;
    end
end

[Nx,Nc] = size(X);
Nk = Nc/3;

%% Sanity checks
if rem(Nk,1); error('3C data required!'); end       % 3-component data only
if isreal(X); X = hilbert(X); end                   % Ensure X is analytic
if nargout > 1 && ~ht && ~h1; h1 = 1; end           % Histograms
Na = min([Na floor(Nx/3)]);                         % Filter order

%% Compute polarization attributes
P.az = zeros(Nx,1,Nk);
P.el = zeros(Nx,1,Nk);
P.in = zeros(Nx,1,Nk);
P.rc = zeros(Nx,1,Nk);
P.pl = zeros(Nx,1,Nk);
W1 = zeros(Nx,1,Nk);
W2 = zeros(Nx,1,Nk);
for k = 1:1:Nk
    x = X(:,3*k-2:3*k);
    C = [x(:,1).*conj(x(:,1)) x(:,1).*conj(x(:,2)) ...
        x(:,1).*conj(x(:,3)) x(:,2).*conj(x(:,1)) ...
        x(:,2).*conj(x(:,2)) x(:,2).*conj(x(:,3)) ...
        x(:,3).*conj(x(:,1)) x(:,3).*conj(x(:,2)) ...
        x(:,3).*conj(x(:,3))];
    
    % This is not exactly averaging, but is equivalent & much faster
    C = circshift(filter((1/Na)*ones(1, Na), 1, C),[-floor(Na/2) 0]);
    [P.pl(:,1,k), P.rc(:,1,k), W2(:,1,k)] = seisjpol(real(C));
    if av
        V1 = zeros(3,Nx);
        for m = 1:1:Nx-Na
            % _______________________________________________________
            % Averaged versions of "instantaneous" parameters
             c = [C(m,1) C(m,2) C(m,3); ...
                  C(m,4) C(m,5) C(m,6); ...
                  C(m,7) C(m,8) C(m,9)];
            % Eigenvector
            [V, d] = eig(c,'nobalance');
            [~, i1]  = sort(diag(d),'descend');
            v1 = V(:,i1(1));
            V1(:,m) = v1;
            W1(m,1,k) = sum(abs(diag(c)));
        end
        [P.az(:,1,k), P.el(:,1,k), P.in(:,1,k)] = seisvpol(V1, z, fca);
    else
        V = [x(:,1)./x(:,3) x(:,2)./X(:,3) ones(Nx,1)]';
        V = V./repmat(sqrt(sum(V.*conj(V))),3,1);
        [P.az(:,1,k), P.el(:,1,k), P.in(:,1,k)] = seisvpol(V, z, fca);
        W1(:,1,k) = sum(real(x).^2,2);
    end
end

% Truncate 
Na = 1+ceil((Na-1)/2);
P.az = P.az(Na:Nx-Na,:,:);
P.el = P.el(Na:Nx-Na,:,:);
P.in = P.in(Na:Nx-Na,:,:);
P.pl = P.pl(Na:Nx-Na,:,:);
P.rc = P.rc(Na:Nx-Na,:,:);
W1 = W1(Na:Nx-Na,:,:);
W2 = W2(Na:Nx-Na,:,:);
P = orderfields(P);

%% Histograms
if ht || h1
    G = gdm(1/dp, repmat(pv,[1 5]), [1 0 1 0 0]);
    
    % Generate weight structure
    W{1} = W1;
    W{2} = W1;
    W{3} = W1;
    W{4} = W2;
    W{5} = W2;

    % Map az, inc to [0, 1]
    if fca; a0 = 2; else a0 = 1; end
    P.az = P.az/((fca+1)*180) + 0.5;
    P.in = P.in/((fca+1)*180) + 0.5;

    if ht
        if h1
            disp(['Warning! h1 == 1 & ht == 1. Only doing ' ...
                ' time slice (ht) histograms.']);
        end
        [H, D, T] = polhist(P, W, 'm', 't', ...
            'G', G, 'L', L, 'La', La, 'tk', 'tk', ...
            'dp', dp, 'hd', hd, 'fca', fca, 'pv', pv);
    else
        [H, D] = polhist(P, W, 'm', '1', ...
            'G', G, 'L', L, 'La', La, 'tk', 'tk', ...
            'dp', dp, 'hd', hd, 'fca', fca, 'pv', pv);
        T = 0;
    end
    P = psqueeze(P);

    % Remap az, in to [-90, 90] or [-180, 180]
    P.az = a0*180*(P.az-0.5);
    P.in = a0*180*(P.in-0.5);
else
    P = psqueeze(P);
end

%% Parse varargout
varargout{1} = P;
if nargout > 1
    varargout{2} = squeeze(W);
    if ht || h1
        varargout{3} = H;
        if nargout > 3
            varargout{4} = D;
            if nargout > 4
                varargout{5} = T;
            end
        end
    end
end

% ================================================================
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

% ================================================================
function [pl, rc, W] = seisjpol(C)
Nx = size(C,1);
W = zeros(Nx, 1);
pl = zeros(Nx, 1);
rc = zeros(Nx, 1);
for m = 1:1:Nx
    [~,u] = eig([C(m,1) C(m,2) C(m,3); ...
                 C(m,4) C(m,5) C(m,6); ...
                 C(m,7) C(m,8) C(m,9)],'nobalance');
    [u,~] = sort(diag(u),'descend');
    rc(m) = 1 - ( (u(2) + u(3)) / (2*u(1)) );
    pl(m) = 1 - ( 2*u(3) / (u(1)+u(2)) );
    
    W(m) = C(1) + C(5) + C(9);
end

% ================================================================
function S = psqueeze(P)
F = fieldnames(P);
for p = 1:1:numel(F)
    P.(F{p}) = squeeze(P.(F{p}));
end
S = P;
