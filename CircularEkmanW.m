%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIRCULAR EKMAN TRANSPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 0.075;
out = calcCircW(-epsilon);
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic;
%%
vlim = [-1.25 -.85];
wlim = [-1 1]*5e-5;

%%
muconts = linspace(-0.1, 0.1, 20);
mvconts = linspace(-.5, -1.25, 20);
wconts = linspace(-3e-5, 3e-5, 20);
wcontsclass = linspace(-1e-5, 1e-5, 40);

subplot(4,2,1);
contourf(X./cr, Y./cr, MUBAR, muconts); shading interp
axis square
colorbar;
set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
grid on
ylabel('M_x');

subplot(4,2, 3)
contourf(X./cr, Y./cr, MVBAR, mvconts); shading interp
axis square
set(gca, 'clim', vlim);
colorbar
set(gca, 'xlim', [-1.25 1.25], 'ylim', [-1.25 1.25]);
grid on
ylabel('M_y');

subplot(4,2,5);
contourf(X./cr, Y./cr, W, wconts); shading interp
axis square
set(gca, 'clim', wlim);
colorbar
set(gca, 'xlim', [-1.25 1.25], 'ylim', [-1.25 1.25]);
grid on
ylabel('W');

subplot(4,2,7);
contourf(X./cr, Y./cr, W-Wclassic, wcontsclass);
axis square
colorbar
set(gca, 'xlim', [-1.25 1.25], 'ylim', [-1.25 1.25]);
grid on
ylabel('W''');

%%
%Cyclone
out = calcCircW(epsilon);
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic;
%%

muconts = linspace(-0.1, 0.1, 20);
mvconts = linspace(-.5, -1.25, 20);
wconts = linspace(-3e-5, 3e-5, 20);
wcontsclass = linspace(-1e-5, 1e-5, 40);

subplot(4,2,2);
contourf(X./cr, Y./cr, MUBAR, muconts); shading interp
axis square
colorbar;
set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
grid on
ylabel('M_x');

subplot(4,2, 4)
contourf(X./cr, Y./cr, MVBAR, mvconts); shading interp
axis square
set(gca, 'clim', vlim);
colorbar
set(gca, 'xlim', [-1.25 1.25], 'ylim', [-1.25 1.25]);
grid on
ylabel('M_y');

subplot(4,2,6);
contourf(X./cr, Y./cr, W, wconts); shading interp
axis square
set(gca, 'clim', wlim);
colorbar
set(gca, 'xlim', [-1.25 1.25], 'ylim', [-1.25 1.25]);
grid on
ylabel('W');

subplot(4,2,8);
contourf(X./cr, Y./cr, W-Wclassic, wcontsclass);
axis square
colorbar
set(gca, 'xlim', [-1.25 1.25], 'ylim', [-1.25 1.25]);
grid on
ylabel('W''');


%%
% wld = wlim(end)/2;
% wconts = linspace(-wld, wld, 20);
% figure
% contourf(X./cr, Y./cr, W-Wclassic, wconts); shading interp
% axis square
% 
% uconts = -(0:.1:.6);
% % uconts = [.5 .5];
% hold on
% [c, h] = contour(X./cr, Y./cr, UMAG, uconts, 'k');
% clabel(c,h)
% hold off
% set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
% set(gca, 'clim', wlim/2);
% grid on
% title(['\epsilon_{sh} = ', num2str(epsilon), '   \epsilon_k = ', num2str(maxom./f)]);