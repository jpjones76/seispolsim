% Figures for test polarization script
figure('Position',[560 200 800 600]);
subplot(4,2,1)
hold on
plot(3 + X(:,1)./c);
plot(2 + X(:,2)./c);
plot(1 + X(:,3)./c);
text(0, 3, 'X$$_z$$', 'interpreter', 'latex', 'horizontalalignment', 'right');
text(0, 2, 'X$$_n$$', 'interpreter', 'latex', 'horizontalalignment', 'right');
text(0, 1, 'X$$_e$$', 'interpreter', 'latex', 'horizontalalignment', 'right');
plot([t1 t1],[0 4],'k--');
plot([t2 t2],[0 4],'k--');
xlim([0 2048]); ylim([0.5 3.5]);
set(gca,'yticklabel','');

subplot(4,2,3)
plot(P.el)
hold on
plot([t1 t1],[0 1],'k--');
plot([t2 t2],[0 1],'k--');
xlim([0 2048]); ylim([0 1]);
ylabel('$$\mathbf\eta$$','interpreter','latex');

subplot(4,2,5)
plot(P.az);
xlim([1 2048]);
hold on
plot([t1 t1],[-pi/2 pi/2],'k--');
plot([t2 t2],[-pi/2 pi/2],'k--');
xlim([0 2048]); ylim([-90 90]);
ylabel('$$\mathbf\theta$$','interpreter','latex');

subplot(4,2,7)
plot(P.in);
hold on
plot([t1 t1],[-pi/2 pi/2],'k--');
plot([t2 t2],[-pi/2 pi/2],'k--');
xlim([0 2048]); ylim([-90 90]);
ylabel('$$\mathbf\phi$$','interpreter','latex');

subplot(4,2,4)
bar(ex, H{2});
set(gca,'xtick',0:0.2:1);
xlim([0 1]);
ylabel('$$\mathbf\eta$$','interpreter','latex');

subplot(4,2,6)
bar(ax, H{1});
set(gca,'xtick',(fca+1)*(-90:30:90));
xlim([-90 90]);
ylabel('$$\mathbf\theta$$','interpreter','latex');

subplot(4,2,8)
bar(ax, H{3});
set(gca,'xtick',(fca+1)*(-90:30:90));
xlim([-90 90]);
ylabel('$$\mathbf\phi$$','interpreter','latex');