function varargout = psiplot(X,S,T,fs,L,varargin)
% h = psiplot(X,S,T,fs,varargin);
% 
% Compute and plot adaptive similarity of polarizations for data sampled
% at fs Hz.
% 
% ====================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2015-12-15

Nk = size(X,2)/3;
tos = 0;
lbfs = 10;
if length(varargin) > 0
    if ~isempty(varargin{1})
        lbs = varargin{1};
    end
    if numel(varargin) > 1
        if ~isempty(varargin{2})
            lbs = varargin{2};
        end
    end
end

% Station and image labels
if ~exist('lbs','var')
    for k = 1:1:Nk
        lbs{k} = num2str(k);
    end
end
lbk = {'$$\mathbf{S_{\theta,k}}$$'; ...
    '$$\mathbf{S_{\eta,k}}$$'; ...
    '$$\mathbf{S_{\phi,k}}$$'; ...
    '$$\mathbf{S_{\nu,k}}$$'; ...
    '$$\mathbf{S_{\rho,k}}$$'};
lbt = {'$$\mathbf{S_{\theta,t}}$$'; ...
    '$$\mathbf{S_{\eta,t}}$$'; ...
    '$$\mathbf{S_{\phi,t}}$$'; ...
    '$$\mathbf{S_{\nu,t}}$$'; ...
    '$$\mathbf{S_{\rho,t}}$$'};

%% Correct times to center
La = T(2)-T(1);
T = (T-La/2)/fs; 
S = circshift(S,[0 round(L/La)]);

% Extract St, Sk
St = S(1:5*Nk,:);
Sk = S(5*Nk+1:10*Nk,:);

% Sk figure
h1 = polsimfig(X, Sk, 'fs', fs, 'td', T, 'pnames', lbk, ...
    'Ls', L/La, 'tos', tos, 'cmp', 1, 'it', 2, ...
    'lbs', lbs, 'lbfs', lbfs);
set(h1(1), 'name', 'Sk relief plot');
drawnow;

% St figure
h2 = polsimfig(X, St, 'fs', fs, 'td', T, 'pnames', lbt, ...
    'Ls', L/La, 'tos', tos, 'cmp', 1, 'it', 2, ...
    'lbs', lbs, 'lbfs', lbfs);
set(h1(1), 'name', 'St relief plot');
drawnow;

% Parse varargout
if nargout
    varargout{1} = [{Sk} {St}];
    if nargout > 1
        varargout{2} = [h1 h2];
    end
end

