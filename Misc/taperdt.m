function Y = taperdt(varargin)
% Y = taperdt(X);
%   Taper and de-mean X. 
%
% Y = taperdt(X,tf);
%   Specify decimal % of X to taper. Default is 0.01.
% 
% Modified from CORAL routines by Ken Creager, kcc@ess.washington.edu
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 0, written before 2010-01-01

tf = 0.005;
X = varargin{1};
if nargin > 1
    tf = varargin{2}/2;
end
[Nx,Nc] = size(X);
m = round(Nx*tf+0.5-eps);
t = (0:2*m-1)'/(2*m-1);
y = .5*(1 - cos(2*pi*t));
tap = [y(1:m); ones(Nx-2*m,1); y(m+1:2*m)];
Y = detrend(X.*repmat(tap,[1 Nc]),'constant');
