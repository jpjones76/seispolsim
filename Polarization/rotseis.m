function Y = rotseis(varargin)
% Y = rotseis(X, Az);
% 
% Rotate matrix of 3 component seismograms X into aziuth vector Az.
%
% X should have column vectors arranged [Z_1 N_1 E_1 Z_2 ... Z_K N_K E_K]
% for a K-sensor array. Treat +Z as down.
%
% For a K-sensor array, if data have Nx samples, X should be size [Nx, K/3]
%
% =====================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 0, written before 2010-01-01

% X and dependent variables
X = varargin{1};
phi = varargin{2};
Nk = size(X,2);
if mod(Nk,3)
    disp('Non-3-component data matrix, exit with error');
    return
else
    Nk = Nk/3;
end

% Repeat if single angle
if numel(phi) == 1
    phi = repmat(phi,[1 Nk]);
end

% Optional inputs.  Reverse rotation if left.
lr = 'r';
if nargin > 2
    lr = varargin{3};
end
if strcmpi(lr,'l')
    phi = - phi;
end

% Rotate
Y = zeros(size(X));
z1 = X(:,1:3:end);
n1 = X(:,2:3:end);
e1 = X(:,3:3:end);

for k = 1:1:Nk
    M = [cosd(phi(k)) sind(phi(k)); ...
        -sind(phi(k)) cosd(phi(k))];
    rt = M*[n1(:,k)'; e1(:,k)'];
    Y(:,3*k-2) = z1(:,k);
    Y(:,3*k-1) = rt(1,:)';
    Y(:,3*k-0) = rt(2,:)';
end
