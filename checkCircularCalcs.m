cl = [-1.25 -.75];

subplot(3,1,1)
pcolor(X, Y, out.MV) ; shading interp
title('Full Analytic Sol');
set(gca, 'clim', cl);
colorbar;

subplot(3,1,2)
pcolor(X, Y, out.MVEPS); shading interp
title('Approx. Analytic Sol');
set(gca, 'clim', cl);
colorbar;

subplot(3,1,3)
pcolor(X, Y, out.MVCLASS); shading interp
title('Classic Sol');
set(gca, 'clim', cl);
colorbar;

%%
figure
subplot(2,1,1)
plot(Y./cr, out.OMEGA(:,100));
hold on
plot(Y./cr,out.ZETA(:,100)); 
plot(Y./cr,out.ZETA(:, 100) - out.OMEGA(:,100));
hold off
grid on
subplot(2,1,2)
plot(Y./cr, (out.MV(:,100))); hold on;
plot(Y./cr, out.MVCLASS(:,100));
plot(Y./cr, out.MVEPS(:,100))
plot(Y./cr, out.MVLIN(:,100)); hold off
% set(gca, 'xlim', [0 200])
grid on

%%
