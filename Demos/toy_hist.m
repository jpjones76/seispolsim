clear all; close all;
setfigdefs;

% Dependent parameters
Nh = 60;
thr = Nh/5;
T = [thr thr].*ones(1,2);
G = gdm(Nh,T,[1 1],{'lin';'gauss'});
lcol = 0.8*ones(1,3);
dcol = 0.2*ones(1,3);
t1 = 0.5:1:Nh-0.5;
dxt = ceil(Nh/4);
xt = [0:dxt:Nh];
yt = [0 5 10];
xl = round(Nh/4);
xu = round(3*Nh/4);
xr = [xl xu];
file1 = 'toy_histograms.eps';
target_dir_win = 'D:\Temp"\';

% Toy histograms
h0 = zeros(Nh,1);

h1a = h0;
h1a(floor(Nh/2)+1) = 5;
h1a(floor(Nh/2)) = 9;
h1b = h0;
h1b = circshift(h1a,[-2 0]);

h0a = circshift(h1a,[-1 0]);
d = hdist([h1a h0a],G{1},'t','chi2');
D0 = d(1,2);

d = hdist([h1a h1b],G{1},'t','chi2');
D1 = d(1,2);

d = hdist([h1a h1b],G{1},'t','qchd');
D2 = d(1,2);

d = hdist([h1a h1b],G{2},'t','qchd');
D3 = d(1,2);

h2a = h1a;
h2b = circshift(h1a,[1-round(0.5*thr) 0]);
d = hdist([h2a h2b],G{2},'t','qchd');
D4 = d(1,2);

% Plot
h0aplot = reshape([zeros(size(h0a))'; h0a'; h0a'; ...
    zeros(size(h0a))'],1,[]);
h1aplot = reshape([zeros(size(h1a))'; h1a'; h1a'; ...
    zeros(size(h1a))'],1,[]);
h1bplot = reshape([zeros(size(h1b))'; h1b'; h1b'; ...
    zeros(size(h1b))' ],1,[]);
h2aplot = reshape([zeros(size(h2a))'; h2a'; h2a'; ...
    zeros(size(h2a))'],1,[]);
h2bplot = reshape([zeros(size(h2b))'; h2b'; h2b'; ...
    zeros(size(h2b))' ],1,[]);
t1plot = reshape([0.5:1:Nh; 0.5:1:Nh; 1.5:1:Nh+1; ...
    1.5:1:Nh+1],1,[]);
fh = figure('PaperPosition',[0.25 0.25 6 8]);

% H0
axes('Position',[0.1 0.82 0.8 0.17]);
hold on;
fill(t1plot,h1aplot,dcol,'Linewidth',2);
fill(t1plot,h0aplot,lcol,'Linewidth',2);
xlim(xr);
ylim([0 1.05*max(yt)]);
text(xu-1, 0.95*max(yt), ...
    ['$${\mathbf D^{(\chi^2)} = ' sprintf('%0.2f',D0) '}$$'], ...
    'fontsize', 14, ...
    'HorizontalAlignment', 'Right', ...
    'VerticalAlignment', 'Cap');
set(gca,'xtick',xt,'xticklabel','','linewidth',2);
set(gca,'ytick',yt);

% H1
axes('Position',[0.1 0.64 0.8 0.17]);
hold on;
fill(t1plot,h1aplot,dcol,'Linewidth',2);
fill(t1plot,h1bplot,lcol,'Linewidth',2);
xlim(xr);
ylim([0 1.05*max(yt)]);
text(xu-1, 0.95*max(yt), ...
    ['$${\mathbf D^{(\chi^2)} = ' sprintf('%0.2f',D1) '}$$'], ...
    'fontsize', 14, ...
    'HorizontalAlignment', 'Right', ...
    'VerticalAlignment', 'Cap');
set(gca,'xtick',xt,'xticklabel','','linewidth',2);
set(gca,'ytick',yt);

% H1 with linear QCHD
axes('Position',[0.1 0.46 0.8 0.17]);
hold on;
plot(t1,max(yt)*circshift(G{1}(1,:),[0 -floor(Nh/2)]),'k--','linewidth',2);
fill(t1plot,h1aplot,dcol,'Linewidth',2);
fill(t1plot,h1bplot,lcol,'Linewidth',2);
xlim(xr);
ylim([0 1.05*max(yt)]);
text(xu-1, 0.95*max(yt), ...
    ['$${\mathbf D^{(QC)} = ' sprintf('%0.2f',D2) '}$$'], ...
    'fontsize', 14, ...
    'HorizontalAlignment', 'Right', ...
    'VerticalAlignment', 'Cap');
set(gca,'xtick',xt,'xticklabel','','linewidth',2);
set(gca,'ytick',yt);

% H1 with Gaussian QCHD
axes('Position',[0.1 0.28 0.8 0.17]);
hold on;
plot(t1,max(yt)*circshift(G{2}(1,:),[0 -floor(Nh/2)]),'k--','linewidth',2);
fill(t1plot,h1aplot,dcol,'Linewidth',2);
fill(t1plot,h1bplot,lcol,'Linewidth',2);
xlim(xr);
ylim([0 1.05*max(yt)]);
text(xu-1, 0.95*max(yt), ...
    ['$${\mathbf D^{(QC)} = ' sprintf('%0.2f',D3) '}$$'], ...
    'fontsize', 14, ...
    'HorizontalAlignment', 'Right', ...
    'VerticalAlignment', 'Cap');
set(gca,'xtick',xt,'xticklabel','','linewidth',2);
set(gca,'ytick',yt);

% H2 with Gaussian QCHD
axes('Position',[0.1 0.1 0.8 0.17]);
hold on;
plot(t1,max(yt)*circshift(G{2}(1,:),[0 -floor(Nh/2)]),'k--','linewidth',2);
fill(t1plot,h2aplot,dcol,'Linewidth',2);
fill(t1plot,h2bplot,lcol,'Linewidth',2);
xlim(xr);
ylim([0 1.05*max(yt)]);
text(xu-1, 0.95*max(yt), ...
    ['$${\mathbf D^{(QC)} = ' sprintf('%0.2f',D4) '}$$'], ...
    'fontsize', 14, ...
    'HorizontalAlignment', 'Right', ...
    'VerticalAlignment', 'Cap');
set(gca,'xtick',xt,'linewidth',2);
set(gca,'ytick',yt);
xlabel('{\bf i [Bin \#]}','FontSize',12);

% Equation gdm2
M = [];
for m = 1:1:Nh;
    M1 = [];
    for n = [1:1:5 Nh-4:1:Nh]
        if n ==5
            M1 = strcat(M1,sprintf(' %4s &','... '));
        else
            if G{1}(m,n)
                M1 = strcat(M1,sprintf(' %4.2f &',G{1}(m,n)));
            else
                M1 = strcat(M1,sprintf(' %4i &',0));
            end
        end
    end
    M1 = M1(1:end-1);
    M1(end+1:end+2) = '\\';
    M = [M; M1];
end

% Print and Copy
print(fh,'-depsc2','-r600','-noui',file1);
if strcmpi(computer,'pcwin64')
    disp(['xcopy ' file1 ' ' target_dir_win ' /Y']);
    system(['xcopy ' file1 ' ' target_dir_win ' /Y']);
end
