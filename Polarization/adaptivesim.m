function S = adaptivesim(X, D, L, La)
% S = adaptivesim(X, D, L, La);
% 
%   Convert distances in cell structure D adaptive similarities using
% STA/LTA ratios of X, where D measures average histogram distance between
% each station k and each station k1 ~= k.
% 
% S = adaptivesim(X, D, L, La);
%
%   Specify t=1 to if using distances computed between adjacent time 
% slice histograms at each station (affects the STA/LTA calculation).
%
% Dependencies: stalta, Dcorr.mat
%
% ======================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.4.1, 2015-12-16

% Required: Data, dist, type
[Nx,Nk] = size(X);
Nk = Nk/3;
Lw  = 4*L;
Dt = D(1:5*Nk,:);
Dk = D(5*Nk+1:10*Nk,:);
St = zeros(size(Dt));
Sk = zeros(size(Dk));
DUt = zeros(Nk, size(Dt,2));
DLt = zeros(Nk, size(Dt,2));

% Get Dmax, Dmin
pstr = {'az';'el';'in';'pl';'rc'};
load('Dcorr.mat','SN','DmaxA','DminA');

%% STA-LTA ratios
if isreal(X) && min(X(:))<0
    X = abs(hilbert(X));
end
Rk = stalta(mean(X,2), L, Lw, La);

Xs = zeros(Nx, Nk);
for k = 1:1:Nk
    Xs(:,k) = sum(X(:,3*k-2:3*k),2);
end
Rt = stalta(Xs, L, Lw, La);

%% Convert D to adaptive similarity
for p = 1:1:5
    D0 = Dt(Nk*(p-1)+1:Nk*p,:);
    D1 = Dk(Nk*(p-1)+1:Nk*p,:);
    Dmax = DmaxA.(pstr{p});
    Dmin = DminA.(pstr{p});

    [~,NS] = size(D);

    % Cross-time similarity
    SN0 = max(-4, 20*log10(Rt(1:size(D,2),:)))';
    for i = 1:1:NS
        for k = 1:1:Nk
            j = max(2, find(SN > SN0(k,i), 1));
            if isempty(j); j = numel(SN)+1; end
            DUt(k,i) = Dmax(j-1);
            DLt(k,i) = Dmin(j-1);
        end
    end
    St(Nk*(p-1)+1:Nk*p,:) = min(1, (1 - 2*min(1, (D0-DLt)./(DUt-DLt))));
    
    % Cross-station similarity
    SN1 = max(-4, 20*log10(Rk(1:size(D,2),:)))';
    DUk = zeros(1, NS);
    DLk = zeros(1, NS);
    for i = 1:1:NS
        j = max(2, find(SN > SN1(i), 1));
        if isempty(j); j = numel(SN)+1; end
        DUk(i) = Dmax(j-1);
        DLk(i) = Dmin(j-1);
    end
    DUk = repmat(DUk, size(D1,1), 1);
    DLk = repmat(DLk, size(D1,1), 1);
    Sk(Nk*(p-1)+1:Nk*p,:) = min(1, (1 - 2*min(1, (D1-DLk)./(DUk-DLk))));
end
S = [St; Sk];
