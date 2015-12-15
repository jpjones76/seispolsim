% Script to generate toy histogram figures
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
xt = [5 10 15 20];
fh = figure('PaperPosition',[0.25 0.25 6 8]);

% H1
axes('Position',[0.1 0.7 0.8 0.19]);
hold on;
fill(t1plot,h1aplot,dcol,'Linewidth',1);
fill(t1plot,h1bplot,lcol,'Linewidth',1);
xlim([0.5 Nh+0.5]);
ylim([0 1.05*max(yt)]);
text(Nh-1, 0.95*max(yt), ...
     ['$${\mathbf D^{(\chi^2)} = ' sprintf('%0.2f',D1) '}$$'], ...
     'HorizontalAlignment', 'Right', ...
     'VerticalAlignment', 'Cap');
set(gca,'xtick',xt,'xticklabel','');
title(['{\bf Histogram Distances, } $$\mathbf{\tau = ' num2str(thr-1) '}$$'],  ...
      'FontSize', 12);
set(gca,'ytick',yt);

% H1 with linear QCHD
axes('Position',[0.1 0.5 0.8 0.19]);
hold on;
plot(t1,max(yt)*circshift(G{1}(1,:),[0 -floor(Nh/2)]),'k--','linewidth',2);
fill(t1plot,h1aplot,dcol,'Linewidth',1);
fill(t1plot,h1bplot,lcol,'Linewidth',1);
xlim([0.5 Nh+0.5]);
ylim([0 1.05*max(yt)]);
text(Nh-1, 0.95*max(yt), ...
     ['$${\mathbf D^{(QC)} = ' sprintf('%0.2f',D2) '}$$'], ...
     'HorizontalAlignment', 'Right', ...
     'VerticalAlignment', 'Cap');
set(gca,'xtick',xt,'xticklabel','');
set(gca,'ytick',yt);

% H1 with Gaussian QCHD
axes('Position',[0.1 0.3 0.8 0.19]);
hold on;
plot(t1,max(yt)*circshift(G{2}(1,:),[0 -floor(Nh/2)]),'k--','linewidth',2);
fill(t1plot,h1aplot,dcol,'Linewidth',1);
fill(t1plot,h1bplot,lcol,'Linewidth',1);
xlim([0.5 Nh+0.5]);
ylim([0 1.05*max(yt)]);
text(Nh-1, 0.95*max(yt), ...
     ['$${\mathbf D^{(QC)} = ' sprintf('%0.2f',D3) '}$$'], ...
     'HorizontalAlignment', 'Right', ...
     'VerticalAlignment', 'Cap');
set(gca,'xtick',xt,'xticklabel','');
set(gca,'ytick',yt);

% H2 with Gaussian QCHD
axes('Position',[0.1 0.1 0.8 0.19]);
hold on;
plot(t1,max(yt)*circshift(G{2}(1,:),[0 -floor(Nh/2)]),'k--','linewidth',2);
fill(t1plot,h2aplot,dcol,'Linewidth',1);
fill(t1plot,h2bplot,lcol,'Linewidth',1);
xlim([0.5 Nh+0.5]);
ylim([0 1.05*max(yt)]);
text(Nh-1, 0.95*max(yt), ...
     ['$${\mathbf D^{(QC)} = ' sprintf('%0.2f',D4) '}$$'], ...
     'HorizontalAlignment', 'Right', ...
     'VerticalAlignment', 'Cap');
set(gca,'xtick',xt);
set(gca,'ytick',yt);
xlabel('{\bf Bin \#}','FontSize',12);
