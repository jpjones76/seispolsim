function varargout = polsimfig(varargin)
% h = polsimfig(X, S, OPTIONS);
%
% Polarization Similarity Time Series
%     Plot a time series of average cross-station polarization similarity
% of X. 
% 
% DEPENDENCIES: setfigdefs.m, AlignYLbl.m
%
% =======================================================================
% Author: Joshua Jones, highly.creative.pseudonym@gmail.com
% Version: 1.0, 2015-12-14

setfigdefs;

% Options
fs      = 100;
cmp     = 3;
xb      = 0;
Np      = 5;
it      = 1;
tos     = 0;
L       = 100;
La      = 50;
ci      = -1; 
ca      = 1;
lbfs    = 8;
colmap = 'jet';
pnames = {  '$$\mathbf{S(\phi)}$$'; ...
            '$$\mathbf{S(\eta)}$$'; ...
            '$$\mathbf{S(\delta)}$$'; ...
            '$$\mathbf{S(\nu)}$$'; ...
            '$$\mathbf{S(\rho)}$$'  };
        
% Required 1: Data
X = varargin{1};
[Nx, Nc] = size(X); 
Nk = Nc/3;
td = 0.5+1/fs : 0.5 : Nx/fs;


% Required 2: Similarity matrix
S = varargin{2};
S(isnan(S)) = 1;

% Optional: Everything else
if nargin > 2
    j = 3;
    while j < nargin
        eval([varargin{j} '= varargin{j+1};']);
        j = j+2;
    end
end

td = td + tos;
t = tos + (1/fs:1/fs:Nx/fs)';
xi = min(td);
xa = max(td);

if xb
    t0 = find(td>=xb, 1); if isempty(t0); t0 = 1; end
    t1 = numel(td)-t0;
    xi = td(t0);
    xa = td(t1);
end

x = X(:,cmp);
xa = min([xa max(t)]);
yy = 1:1:Nk;

switch lower(cmp)
    case 1
        cstr = 'z';
    case 2
        cstr = 'n';
    otherwise
        cstr = 'e';
end

% Axes positions
yh = 0.14;
yah = yh - 0.01;
pos = zeros(6,4);
wid = 0.73;
xl = 0.12;
if it == 1
    wid = wid + 0.08;
    xl = xl + 0.03;
end
for k = 1:1:6
    pos(k,:) = [xl 0.1+(k-1)*yh wid yah];
end
pos(:,2) = flipud(pos(:,2));

% Create figure
h = zeros(1,7);
h(1) = figure('PaperPosition',[0.5 0.5 8 10]);

% ____________________________________________________
% Seismogram
h(2) = axes('Position',pos(1,:));
plot(t, x, 'k-', 'linewidth', 1);
ylabel(['$$\mathbf{x_{1' cstr '}}$$'], ...
    'interpreter','latex',...
    'FontSize',11, ...
    'rotation',0, ...
    'horizontalalignment','right', ...
    'verticalalignment','middle');

% Axis manipulation
ym = 1.1*max(abs(X(:,1)));
set(gca, ...
    'tickdir','out', ...
    'xticklabel',{}, ...
    'xlim', [xi xa], ...
    'ylim', [-ym ym]);


% ____________________________________________________
% Images for polarization similarity
switch it
    case 1       
        for p = 1:5
            
            % Create axes
            h(p+2) = axes('Position',pos(p+1,:));          
            
            % Similarities for plot
            s = S(p,:);
            
            % Create image
            plot(td, s,'k-','linewidth',2);
            yi = floor(10*min(s(:,t0:t1)))/10; if yi == 1; yi = 0.9; end
            ya = ceil(10*max(s(:,t0:t1)))/10;
            
            % Label
            ylabel(pnames{p}, ...
                'FontSize',11, ...
                'rotation',0, ...
                'horizontalalignment','right', ...
                'verticalalignment','middle');
            
            % Axis manipulation
            dy = 0.1;
            if ya - yi >= 1; dy = 0.4; elseif ya - yi >= 0.5; dy = 0.2; end
            set(h(p+2), ...
                'fontsize', 11, ...
                'position', pos(p+1,:), ...
                'xlim',[xi xa], ...
                'ylim',[yi ya], ...
                'ytick',fliplr(ya:-dy:yi), ...
                'tickdir', 'out');
            if p < Np
                set(gca,'xticklabel', {});
            end
            if max(max(s)) > 0.99 && ya == 1
                set(gca,'ylim',[yi 1.05]);
            end
        end
        AlignYLbl(h(1), 0.12, 'center');
        
    case 2
        axes(h(2));
        text(xi+1.02*(xa-xi), 0, ['$$\mathbf{v\ \bigl [ \frac{m}{s} \bigr ]}$$'], ...
            'FontSize', 11, ...
            'rotation', 0, ...
            'horizontalalignment','left', ...
            'verticalalignment','middle');
        
        for p = 1:5
            % Create axes
            h(p+2) = axes('Position',pos(p+1,:));
            
            % Similarities for imagesc
            s = S(1+(p-1)*Nk:p*Nk,:);
            
            % Moving average filter
            R = L/La;
            if ~mod(R,2)
                R = R+1;
            end
            sfilt = ones(1,floor(R))/R;
            s2 = circshift(filter(sfilt, 1, s'),[0 -floor(R/2)])';
            clims = [ci ca];
            
            % Create image
            imagesc(td, yy, s, clims);
            colormap(colmap);
            
            % Color bar
            cbh(p) = colorbar;
            
            set(cbh(p), ...
                'Position',[0.01+pos(p+1,1)+pos(p+1,3) pos(p+1,2) 0.025 yah], ...
                'ytick',-1:0.5:1, ...
                'ylim', [ci-eps ca+eps], ...
                'yaxislocation', 'right', ...
                'FontSize', 9);
            
            % Label colorbar
            ylabel(cbh(p), pnames{p}, ...
                'interpreter','latex',...
                'FontSize',11, ...
                'rotation',0, ...
                'horizontalalignment','left', ...
                'verticalalignment','middle');
            
            % Reset axes position, axis manipulation
            set(h(p+2), ...
                'position', pos(p+1,:), ...
                'ytick', yy, ...
                'yticklabel', '', ...
                'xlim',[xi xa], ...
                'tickdir', 'out');
            
            % Text labels
            axes(h(p+2));
            for k = 1:1:Nk
                text(-0.01*(xa-xi)+xi, k, lbs(k), ...
                    'interpreter','latex',...
                    'FontSize', lbfs, ...
                    'rotation', 0, ...
                    'horizontalalignment','right', ...
                    'verticalalignment','middle');
            end
            
            if p < Np
                set(gca,'xticklabel', {});
            end
        end
        
        AlignYLbl(h(1), 0.08, 'center');
    otherwise
end



axes(h(end));
xlabel('{\bf Time [s]}','interpreter','latex');
drawnow;

if nargout
    varargout{1} = h;
end
