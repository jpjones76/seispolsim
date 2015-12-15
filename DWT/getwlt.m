function varargout = getwlt(wlt, J, mo)
% [wf, Lj, Sh] = getwlt(wlt, J, mo);
% [wf, Lj, Sh, Seq] = getwlt(wlt, J, mo);
%      Load wavelet filter; optionally, compute Lj, Sh, Seq by the center
% of energy method. Formulas are from Wickerhauser (1994), Hess-Nielsen and
% Wickerhauser (1996), and Percival & Walden (2000).
% 
% wf    Wavelet filter. 1x2 cell structure. 
%           wf{1} = scaling (low-pass) filter g. 
%           wf{2} = wavelet (high-pass) filter h.
% Lj    Number of coefficients affected by circular boundary condition. 
% Sh    Circular shifts to align W in time. 
% Seq   Sequency vectors that determine whether one filters with g or h.
%       if Seq{j,n}(end) = 0, filter with g; else, use h.
% 
% NOTES
% 1)  Lj(j)+1 is the first coefficient unaffected by the circular boundary
% condition at wavelet level j.
% 
% 2)  All shifts are negative because they're advances; to align W{j,n},
% use circshift(W{j,n},[Sh(j,n) 0]), not circshift(W{j,n},[-Sh(j,n) 0]).
%
% DEPENDENCIES: wscaledb.mat
%
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.1, Last modified: 2015-12-14

%% Get wavelet filter
load wscaledb.mat;
eval(['g = wscale.' wlt ';']);         % Load scaling filter 
L = length(g);
if mo
    g = g./sqrt(2);                     % Convert to "maximal overlap"
end
h = zeros(size(g));
for l = 1:L
    h(l) = g(L-l+1)*((-1)^(l-1));       % Create wavelet filter 
end
wf = [{g} {h}];
varargout{1} = wf;                      % Send to varargout

%% Calculate effective filter width Lj
if nargout > 1
    Lj = zeros(1, J);
    for j = 1:J
        if mo
            Lj(j) = ((2^j)-1)*(L-1);            % P&W (198a)
        else
            Lj(j) = ceil((L-2)*(1-(1/2^j)));    % P&W (146a)
        end
    end
    varargout{2} = Lj;
end

%% Calculate shifts to make wavelet coeffs zero phase
if nargout > 2
    % Initialize sequency cell structure
    seq = cell(J,2^J);
    sh = zeros(J, 2^J);
    
    % Compute center of energy for each filter
    eg = dot(wf{1}.^2,(0:1:L-1));
    eh = dot(wf{2}.^2,(0:1:L-1));
    
    for j = 1:J
        for n = 1:2^j
            % Update sequency; indexing recalculated from P&W section 6.2
            if mod(n,4) < 2
                seq{j,n} = [seq{j,n} 0];
            else
                seq{j,n} = [seq{j,n} 1];
            end
            
            % Compute shifts from sequency
            sh1 = dot(seq{j,n},(2.^(0:1:j-1)));
            sh0 = dot(1-seq{j,n},(2.^(0:1:j-1)));
            sh(j,n) = sh0*eg + sh1*eh;
            if j < J
                seq{j+1,(2*n)-1} = seq{j,n};
                seq{j+1,2*n} = seq{j,n};
            end
        end
    end
    varargout{3} = -1.*round(sh);
    if nargout > 3
        varargout{4} = seq;
    end
end
