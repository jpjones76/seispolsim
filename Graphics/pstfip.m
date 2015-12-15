function ah = pstfip(X, Wd, S, T, fs)
% h = pstfip(X, Wd, S, T, fs)
% 
%       Pol. Sim. time-frequency intensity plot using MODWT detail
%  coefficients Wd, with polarization similarity S, computed from data X
%  sampled at frequency fs, with each histogram i starting at sample
%  numbers T(i).
% 
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, last modified 2015-12-15

Np = 5;
[Nx,Nk] = size(X);
Nk = Nk/3;
cmp = 1;
cstr = 'z';
J = size(Wd,1);
T = T/fs;

pnames{1} = {'$$\mathbf{S_{k,\theta}}$$'; ...
    '$$\mathbf{S_{k,\eta}}$$'; ...
    '$$\mathbf{S_{k,\phi}}$$'; ...
    '$$\mathbf{S_{k,\nu}}$$'; ...
    '$$\mathbf{S_{k,\rho}}$$'};
pnames{2} = {'$$\mathbf{S_{t,\theta}}$$'; ...
    '$$\mathbf{S_{t,\eta}}$$'; ...
    '$$\mathbf{S_{t,\phi}}$$'; ...
    '$$\mathbf{S_{t,\nu}}$$'; ...
    '$$\mathbf{S_{t,\rho}}$$'};
pnames{3} = {'$$\mathbf{\bar{S}_{k,\theta}}$$'; ...
    '$$\mathbf{\bar{S}_{k,\eta}}$$'; ...
    '$$\mathbf{\bar{S}_{k,\phi}}$$'; ...
    '$$\mathbf{\bar{S}_{k,\nu}}$$'; ...
    '$$\mathbf{\bar{S}_{k,\rho}}$$'};
pnames{4} = {'$$\mathbf{\bar{S}_{t,\theta}}$$'; ...
    '$$\mathbf{\bar{S}_{t,\eta}}$$'; ...
    '$$\mathbf{\bar{S}_{t,\phi}}$$'; ...
    '$$\mathbf{\bar{S}_{t,\nu}}$$'; ...
    '$$\mathbf{\bar{S}_{t,\rho}}$$'};

f{1} = 1/Nk:1/Nk:J;
f{2} = f{1};
f{3} = 1:1:J;
f{4} = f{3};
yh{1} = 0.5/Nk;
yh{2} = yh{1};
yh{3} = 0.5;
yh{4} = yh{3};

ft = cell(1,4);
fl = cell(1,4);
for i = 1:1:4
    ft{i} = 0:2:J;
    for n = 1:1:numel(ft{i})
        fl{i}{n} = sprintf('%0.2f',fs/2^(J-ft{i}(n)+1));
    end
    ft{i} = ft{i}+yh{i};
end

% Fill M
for j = 1:1:J
    n = 2;
    if ~isempty(S{j,n})
        s = S{j,n};
        
        er0 = (J-j + mod(n-1,2));
        sr0 = er0;
        er = Nk*er0;
        sr = er - Nk + 1;
        nrep = 0;
        for p = 1:1:Np
            St = s(1+(p-1)*Nk:p*Nk,:);
            Sk = s(1+Nk*(Np+p-1):(Np+p)*Nk,:);
            if nrep
                M{p,1}(sr:er,:) = repmat(Sk,[nrep 1]);
                M{p,2}(sr:er,:) = repmat(St,[nrep 1]);
                M{p,3}(sr0:er0,:) = repmat(mean(Sk),[er0-sr0+1 1]);
                M{p,4}(sr0:er0,:) = repmat(mean(St),[er0-sr0+1 1]);
            else
                M{p,1}(sr:er,:) = Sk;
                M{p,2}(sr:er,:) = St;
                M{p,3}(sr0:er0,:) = mean(Sk);
                M{p,4}(sr0:er0,:) = mean(St);
            end
        end
    end
end

% ______________________________________________
% Create figures

% Axes positions
yht = 0.15;
yah = yht-0.01;
pos = zeros(6,4);
for i = 1:1:6
    pos(i,:) = [0.15 0.1+(i-1)*yht 0.7 yah];
end
pos(:,2) = flipud(pos(:,2));

% Frequencies & Names
figname{1} = 'Sk';
figname{2} = 'St';
figname{3} = 'Sk avg';
figname{4} = 'St avg';
ah = zeros(4,7);

for g = 1:1:4
    ah(g,1) = figure('PaperPosition',[0.5 0.5 8 10], ...
        'Name', figname{g});
    
    % ______________________________________________
    % Seismogram
    ah(g,2) = axes('Position',pos(1,:));
    tx = 1/fs:1/fs:Nx/fs;
    plot(tx,X(:,cmp),'k-');
    
    ylabel('$$\mathbf{v}$$', ...
        'interpreter', 'latex', ...
        'FontSize',11, ...
        'rotation',0, ...
        'horizontalalignment','right', ...
        'verticalalignment','middle');
    
    % Axis manipulation
    ym = 1.1*max(abs(X(:,1)));
    set(gca, ...
        'tickdir','out', ...
        'xticklabel',{}, ...
        'xlim', [min(T) max(T)], ...
        'ylim', [-ym ym]);
    text(1.02*max(T), 0, ['$$\mathbf{x_{1,' cstr '}(t)}$$'], ...
        'interpreter', 'latex', ...
        'FontSize',11, ...
        'rotation',0, ...
        'horizontalalignment','left', ...
        'verticalalignment','middle');
    
    % ______________________________________________
    % Contour plots: distances between histograms
    for p = 1:5
        % Create axes
        ah(g,p+2) = axes('Position',pos(p+1,:));
        
        % Similarities for imagesc
        sim = M{p,g};
        imagesc(T,f{g},sim,[-1 1]);
        colormap('jet');
        axis xy;

        % Label
        ylabel('\bf{F}', ...
            'interpreter', 'latex', ...
            'FontSize',10, ...
            'rotation',0, ...
            'horizontalalignment','right', ...
            'verticalalignment','middle');
        
        % Color bar
        cbh = colorbar;
        set(cbh,'Position',[0.01+pos(p+1,1)+pos(p+1,3) pos(p+1,2) 0.025 yah], ...
            'ytick', [-1:0.4:1], ...
            'ylim', [-1-eps 1+eps], ...
            'yaxislocation', 'right', ...
            'FontSize',9);
        ylabel(cbh, pnames{g}{p}, ...
            'interpreter', 'latex', ...
            'FontSize',11, ...
            'rotation',0, ...
            'horizontalalignment','left', ...
            'verticalalignment','middle');
        
        % Reset axes position, axis manipulation
        set(ah(g,p+2), ...
            'position', pos(p+1,:), ...
            'ytick', ft{g}, ...
            'yticklabel', fl{g}, ...
            'xlim',[min(T) max(T)], ...
            'ylim',[yh{g} max(f{g})+yh{g}], ...
            'tickdir', 'out');
        if p < Np
            set(gca,'xticklabel', {});
        end
    end
    xlabel('{\bf Time [s]}');
    drawnow;
end
