%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIRCULAR EKMAN TRANSPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 0.075;
cr = 60e3;
out = calcCircW(-epsilon);
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic; UMAG= out.UMAG;
%%
% vlim = [-1.25 -.85];
% wlim = [-1 1]*51e-5;
% wlim = [-1 1].*1e-5;


%%
figure
nc = 24;16;
muconts = linspace(-0.2, 0.2, nc);
uclim = muconts([1 end]);
muconts = [-1 muconts 1];
mvconts = linspace(-.5, -1.5, nc);
% mvconts = linspace(-.75, -1.5, nc);
vclim = fliplr(mvconts([1 end]));
mvconts = [-2 mvconts 0];
wconts = linspace(-5e-5, 5e-5, nc);
% wconts = linspace(-1e-5, 1e-5, nc);
wclim = wconts([1 end]);
wconts = [-1 wconts 1];
wcontsclass = linspace(-2e-5, 2e-5, nc);
uconts = -1:.1:1;
% uconts = -[max(max(UBAR)) max(max(UBAR))];
xl = [-1.25 1.25];
xl = [-1.1 1.1];

gap = [.05 .01]; margh = .1; margw = .1;

subtightplot(4,2,1, gap, margh, margw);
contourf(X./cr, Y./cr, MUBAR, muconts); shading interp
axis square
colorbar;
set(gca, 'xlim', xl, 'ylim', xl);
% grid on
ylabel('M_x');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
plot(cos(0:.1:2*pi), sin(0:.1:2*pi), 'k'); 
plot(0.5.*cos(0:.1:2*pi), 0.5.*sin(0:.1:2*pi), 'k'); 
plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off

set(gca, 'clim', uclim);
title('Anticyclone')

subtightplot(4,2,3, gap, margh, margw);
contourf(X./cr, Y./cr, MVBAR, mvconts); shading interp
axis square
% set(gca, 'clim', vlim);
colorbar
set(gca, 'xlim', xl, 'ylim', xl);
% grid on
ylabel('M_y');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
plot(cos(0:.1:2*pi), sin(0:.1:2*pi), 'k'); 
plot(0.5.*cos(0:.1:2*pi), 0.5.*sin(0:.1:2*pi), 'k'); 
plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
set(gca, 'clim',vclim );

subtightplot(4,2,5, gap, margh, margw);
contourf(X./cr, Y./cr, W, wconts); shading interp
axis square
% set(gca, 'clim', wlim);
colorbar
set(gca, 'xlim', xl, 'ylim', xl);
grid on
ylabel('W');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
plot(cos(0:.1:2*pi), sin(0:.1:2*pi), 'k'); 
plot(0.5.*cos(0:.1:2*pi), 0.5.*sin(0:.1:2*pi), 'k'); 
plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
set(gca, 'clim', wclim);

subtightplot(4,2,7, gap, margh, margw);
contourf(X./cr, Y./cr, W-Wclassic, wcontsclass);
axis square
colorbar
set(gca, 'xlim', xl, 'ylim', xl);
grid on
ylabel('W''');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
plot(cos(0:.1:2*pi), sin(0:.1:2*pi), 'k'); 
plot(0.5.*cos(0:.1:2*pi), 0.5.*sin(0:.1:2*pi), 'k'); 
plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
set(gca, 'clim', wcontsclass([1 end]));

%%
%Cyclone
out = calcCircW(epsilon);
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic;
%%

subtightplot(4,2,2, gap, margh, margw);
contourf(X./cr, Y./cr, MUBAR, muconts); shading interp
axis square
colorbar;
set(gca, 'xlim', xl, 'ylim', xl);
% grid on
ylabel('M_x');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
plot(cos(0:.1:2*pi), sin(0:.1:2*pi), 'k'); 
plot(0.5.*cos(0:.1:2*pi), 0.5.*sin(0:.1:2*pi), 'k'); 
plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
set(gca, 'clim', uclim);
title('Cyclone')

subtightplot(4,2,4, gap, margh, margw);
contourf(X./cr, Y./cr, MVBAR, mvconts); shading interp
axis square
% set(gca, 'clim', vlim);
colorbar
set(gca, 'xlim', xl, 'ylim', xl);
% grid on
ylabel('M_y');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
plot(cos(0:.1:2*pi), sin(0:.1:2*pi), 'k'); 
plot(0.5.*cos(0:.1:2*pi), 0.5.*sin(0:.1:2*pi), 'k'); 
plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
set(gca, 'clim',vclim );

subtightplot(4,2,6, gap, margh, margw);
contourf(X./cr, Y./cr, W, wconts); shading interp
axis square
% set(gca, 'clim', wlim);
colorbar
set(gca, 'xlim', xl, 'ylim', xl);
grid on
ylabel('W');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
plot(cos(0:.1:2*pi), sin(0:.1:2*pi), 'k'); 
plot(0.5.*cos(0:.1:2*pi), 0.5.*sin(0:.1:2*pi), 'k'); 
plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
set(gca, 'clim', wclim);

subtightplot(4,2,8, gap, margh, margw);
contourf(X./cr, Y./cr, W-Wclassic, wcontsclass);
axis square
colorbar
set(gca, 'xlim', xl, 'ylim', xl);
grid on
ylabel('W''');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
plot(cos(0:.1:2*pi), sin(0:.1:2*pi), 'k'); 
plot(0.5.*cos(0:.1:2*pi), 0.5.*sin(0:.1:2*pi), 'k'); 
plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
set(gca, 'clim', wcontsclass([1 end]));

colormap(flipud(othercolor('RdBu11')));
set(gcf,'Color', 'w', 'Position', [389   212   727   764])

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