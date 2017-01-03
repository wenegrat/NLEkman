%%
% Compare theoretical solutions
%%
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
f = 2*2*pi./86400*sind(30); 1e-4;
rho = 1035;

%%
X = ncread(statefile, 'X');
Y = ncread(statefile, 'Y');
Z = ncread(statefile, 'Z');
Zl = ncread(statefile, 'Zl');
T = ncread(diagfile, 'T');

nx = length(X); ny = length(Y); nz = length(Z);
dx = X(2)-X(1)
dy = Y(2)-Y(1)
dz = Z(1)-Z(2) %surface only
ts = T(2)-T(1)
nt = length(T);
slice = {0, 0, 0, [nt-100 nt]};

%%
U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, slice);
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, slice);
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, slice);
TAU = GetVar(statefile, etanfile, {'oceTAUX', '(1)'}, {slice{1},slice{2}, [1 1], slice{4}});
 ETA = GetVar(statefile, etanfile, {'ETAN', '(1)'}, {0, 0, [1 1], 0});
%%
navg = 20;

Uavg(:,:,:,1) = nanmean(U(:,:,:,end-navg:end),4);
Vavg(:,:,:,1) = nanmean(V(:,:,:,end-navg:end),4);

Tavg = nanmean(TAU(:,:,:,end-navg:end),4);

ETAavg = nanmean(ETA(:,:,:,end-navg:end), 4);
%%
[omega, phi,theta, UB, VB] = calculateOmega(Uavg, Vavg, dx, dy, f);

% ZETAF = DPeriodic(V, dx, 'x') - DPeriodic(U, dy, 'y');
% TEMP
% TAU = .1.*ones(size(ZETA));
Up = Uavg - repmat(Uavg(:,:,end,:), [1 1 nz 1]);
Vp = Vavg - repmat(Vavg(:,:,end,:), [1 1 nz 1]);

Upt = sum(Up, 3).*dz;
Vpt = sum(Vp, 3).*dz;
Upx = DPeriodic(Upt, dx, 'x');
Vpy = DPeriodic(Vpt, dy, 'y');
WDIA = Upx + Vpy;
%  WDIA2 = repmat(nanmean(WDIA(:,:,:,end-navg:end),4), [1 1 1 1]);
%  WDIA(:,:,1,end) = WDIA2;

ZETA = DPeriodic(Vavg(:,:,end,:), dx, 'x') - DPeriodic(Uavg(:,:,end,:), dy, 'y');
% STERN SOLUTIONS

VTRANS_STERN = -Tavg./(rho*(f+ZETA));
% VTRANS_STERN = -TAU./(rho*(f+ZETA));

WSTERN = DPeriodic(VTRANS_STERN, dy, 'y');

% NIILER SOLUTIONS
[WNI] =calculateNISols(Uavg, Vavg, dx, f,squeeze(Tavg(1,1)));


%% CURVED SOLUTIONS
% theta = phi-pi/2;

denom = ( (f + 2.*omega).*(f+ZETA) - omega.^2);
% denom = f^2 +f.*ZETA + f.*2*omega + 2*omega.*ZETA - omega.^2;

VTRANS = - (f + omega + 2.*omega.*sin(theta).^2 + ZETA.*cos(theta).^2)./denom.*Tavg./rho;
UTRANS =  (ZETA-2.*omega)./denom.*Tavg./rho.*sin(theta).*cos(theta);

%Approximate versions
 UTRANSA = (ZETA-2.*omega).*sin(theta).*cos(theta).*Tavg./(f.^2.*rho);
%  UTRANSA = (ZETA-2.*omega).*(1 - ZETA-2*omega).*sin(theta).*cos(theta).*TAU./(f.^2.*rho);
 VTRANSA = - (f - (ZETA-omega).*sin(theta).^2 - omega.*cos(theta).^2).*Tavg./(f.^2.*rho);

%  UTRANS(97:103,92:98,:,:) = NaN;
%  VTRANS(97:103,92:98,:,:) = NaN;
[Vy, Vx] = gradient(VTRANS, dx);
[Uy, Ux] = gradient(UTRANS, dx);

WCURV =  Ux+Vy;
% WCURV = Vy;
% WCURV( 95:105,90:100,:,:) = NaN;
%% MAIN PLOT
xl = [300 700];
yl = [275 675];
% tim =480;
L = 75e3;
epsilon = 0.158/(f*L);
fs = 20;
fs2 = 12;
% max(max(squeeze(UB(:,:,1,tim))))./(f*L);
tau = .1./1035;

gap = [.1 .01]; margh=.1; margw=.1;
nc = 30;
figure
subtightplot(6,3, [1 4], gap, margh, margw)
cl = linspace(-1, 1,nc);
contourf(X./1000, Y./1000, squeeze(Upt(:,:)).'/(epsilon*tau./f), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{M_{x}^{NUM}}{\epsilon\tau_o/(\rho f)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(305, 650, 'A)', 'FontSize', fs2);

subtightplot(6,3, [7 10], gap, margh, margw)
cl = linspace(-1.05, -.95,nc);
contourf(X./1000, Y./1000, squeeze(Vpt(:,:)).'/(tau./f), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{M_{y}^{NUM}}{\tau_o/(\rho f)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(305, 650, 'B)', 'FontSize', fs2);

subtightplot(6,3,[13 16], gap, margh, margw)
cl = linspace(-0.2, 0.2,nc);
cl = linspace(-6, 6,nc);

contourf(X./1000, Y./1000, squeeze(WDIA(:,:)).'/(epsilon.*tau./(f*L)), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{w_{e}^{NUM}}{\epsilon\tau_o/(\rho f R)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(305, 650, 'C)', 'FontSize', fs2);

subtightplot(6,3, [2 3 5 6 8 9], gap, margh, margw)
% cl = linspace(-0.2, 0.2,nc);
contourf(X./1000, Y./1000, squeeze(WDIA(:,:)-WSTERN(:,:)).'/(epsilon.^2.*tau./(L*f)), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{w_e^{NUM} - w_e^{STERN}}{\epsilon^2\tau_o/(\rho f R)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(305, 655, 'D)', 'FontSize', fs);

subtightplot(6,3, [11 12 14 15 17 18], gap, margh, margw)
% cl = linspace(-0.25, 0.25,nc);
contourf(X./1000, Y./1000, squeeze(WDIA(:,:)-WCURV(:,:)).'/(epsilon.^2.*tau./(L*f)), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
hold on
d = 30;
rectangle('Position', [496-d/2 480-d/2 d d], 'Curvature', [1 1], 'FaceColor', 'w');
hold off
colormap(flipud(othercolor('RdBu11')));
set(gca, 'FontSize', fs2);
t = text(305, 655, 'E)', 'FontSize', fs);
% a=annotation('textbox',[0.5 0.5 1 1],'String','A)','FitBoxToText','on');
title('$\frac{w_e^{NUM} - w_e^{CURVE}}{\epsilon^2\tau_o/(\rho f R)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)

set(gcf, 'Color', 'w', 'Position', [267          22        1135         952]);

%% MAIN PLOT 6 PANEL
xl = [300 700];
yl = [275 675];
% tim =480;
L = 75e3;
epsilon = 0.158/(f*L);
fs = 20;
fs2 = 12;
% max(max(squeeze(UB(:,:,1,tim))))./(f*L);
tau = .1./1035;

gap = [.2 .01]; margh=.1; margw=.1;
nc = 30;
figure

subtightplot(2,3, 1, gap, margh, margw)
cl = linspace(-1, 1,nc);
contourf(X./1000, Y./1000, squeeze(ETAavg(:,:)).'/(max(max(ETAavg))), cl);
% mg = sqrt(UB.^2 + VB.^2);
% contourf(X./1000, Y./1000, squeeze(mg).'/(max(max(abs(mg)))), cl);
contourf(X./1000, Y./1000, squeeze(ZETA(:,:)).'/(max(max(abs(ZETA)))), cl);

grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{\zeta}{f}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(315, 635, 'A)', 'FontSize', fs2);
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');

subtightplot(2,3, 4, gap, margh, margw)
cl = linspace(-1, 1,nc);
contourf(X./1000, Y./1000, squeeze(Upt(:,:)).'/(epsilon*tau./f), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{M_{x}^{NUM}}{\epsilon\tau_o/(\rho f)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(315, 635, 'B)', 'FontSize', fs2);
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');

subtightplot(2,3, 2, gap, margh, margw)
cl = linspace(-1.05, -.95,nc);
contourf(X./1000, Y./1000, squeeze(Vpt(:,:)).'/(tau./f), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{M_{y}^{NUM}}{\tau_o/(\rho f)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(315, 635, 'C)', 'FontSize', fs2);
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');

subtightplot(2,3,5, gap, margh, margw)
cl = linspace(-0.2, 0.2,nc);
cl = linspace(-6, 6,nc);

contourf(X./1000, Y./1000, squeeze(WDIA(:,:)).'/(epsilon.*tau./(f*L)), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{w_{e}^{NUM}}{\epsilon\tau_o/(\rho f R)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(315, 635, 'D)', 'FontSize', fs2);
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');

subtightplot(2,3, 3, gap, margh, margw)
% cl = linspace(-0.2, 0.2,nc);
contourf(X./1000, Y./1000, squeeze(WDIA(:,:)-WSTERN(:,:)).'/(epsilon.^2.*tau./(L*f)), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
set(gca, 'FontSize', fs2);
title('$\frac{w_e^{NUM} - w_e^{STERN}}{\epsilon^2\tau_o/(\rho f R)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)
t = text(315, 635, 'E)', 'FontSize', fs2);
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');

subtightplot(2,3,6, gap, margh, margw)
% cl = linspace(-0.25, 0.25,nc);
contourf(X./1000, Y./1000, squeeze(WDIA(:,:)-WCURV(:,:)).'/(epsilon.^2.*tau./(L*f)), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
xlabel('x (km)', 'Interpreter', 'latex', 'FontSize', fs2);
ylabel('y (km)', 'Interpreter', 'latex', 'FontSize', fs2);
hold on
d = 30;
rectangle('Position', [496-d/2 480-d/2 d d], 'Curvature', [1 1], 'FaceColor', 'w');
hold off
colormap(flipud(othercolor('RdBu11')));
set(gca, 'FontSize', fs2);
t = text(315, 635, 'F)', 'FontSize', fs2);
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');

% a=annotation('textbox',[0.5 0.5 1 1],'String','A)','FitBoxToText','on');
title('$\frac{w_e^{NUM} - w_e^{CURVE}}{\epsilon^2\tau_o/(\rho f R)}$', 'Interpreter', 'Latex', 'FontSize', fs, 'Rotation', 0)

set(gcf, 'Color', 'w', 'Position', [ 303         367        1124         550]);

%% PLOT ERRORS IN APPROXIMATION
figure
subplot(2,2,1)
cl = linspace(-1, 1,nc)*5;
contourf(X./1000, Y./1000, squeeze(Upt(:,:)-UTRANSA(:,:)).'/(epsilon.^2*tau./f), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
subplot(2,2,2)
cl = linspace(-1, 1,nc)*5;
contourf(X./1000, Y./1000, squeeze(UTRANS(:,:)).'/(epsilon.^2*tau./f), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square

subplot(2,2,3)
cl = linspace(-1, 1,nc)*5;
contourf(X./1000, Y./1000, squeeze(Vpt(:,:)-VTRANSA(:,:)).'/(epsilon.^2*tau./f), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
subplot(2,2,4)
cl = linspace(-1, 1,nc)*5;
contourf(X./1000, Y./1000, squeeze(Vpt(:,:)-VTRANS_STERN(:,:)).'/(epsilon.^2*tau./f), cl);
grid on
colorbar
set(gca, 'clim', cl([1 end]));
set(gca, 'xlim', xl, 'ylim', yl);
axis square
%%
mask = ones(size(WCURV));
mask(92:108, 88:96,:,:) = NaN;
werr_stern = squeeze(WDIA(:,:,1,end) - WSTERN(:,:,1,end)).*mask(:,:,1,end);
werr_ap = squeeze(WDIA(:,:,1,end) - WCURV(:,:,1,end)).*mask(:,:,1,end);

wesv = reshape(werr_stern, nx*ny, 1);
weav = reshape(werr_ap, nx*ny, 1);
%  UTRANS(98:102,93:97,:,:) = NaN;
%  VTRANS(98:102,93:97,:,:) = NaN;

% nanstd(wesv)
% nanstd(weav)

sqrt(nanmean(wesv.^2))
sqrt(nanmean(weav.^2))

max(abs(wesv))
max(abs(weav))


%%

cl = [-1 1].*1e-6;
cle = cl./10;

dep = 1;
tim = 480;
figure
subplot(4,2,1)
pcolor(X./1000, Y./1000, squeeze(WDIA(:,:,dep, tim)).');
shading interp
colorbar
set(gca, 'clim', cl);
title('W_{NUMERIC}');
% subplot(3,2,2)
% pcolor(X./1000, Y./1000, squeeze(WCURV(:,:,1, tim) - WSTERN(:,:,1, tim)).');
% shading interp
% colorbar
% set(gca, 'clim', cle);

subplot(4,2,3)
pcolor(X./1000, Y./1000, squeeze(WSTERN(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cl);
title('W_{STERN}');

subplot(4,2,4)
pcolor(X./1000, Y./1000, squeeze(WDIA(:,:,dep, tim) - WSTERN(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cle);
title('W_{NUMERIC}-W_{STERN}');

subplot(4,2,5)
pcolor(X./1000, Y./1000, squeeze(WNI(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cl);
title('W_{NI}');

subplot(4,2,6)
pcolor(X./1000, Y./1000, squeeze(WDIA(:,:,dep, tim) - WNI(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cle);
title('W_{NUMERIC}-W_{NI}');

subplot(4,2,7)
pcolor(X./1000, Y./1000, squeeze(WCURV(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cl);
title('W_{CURV}');

subplot(4,2,8)
pcolor(X./1000, Y./1000, squeeze(WDIA(:,:,dep, tim) - WCURV(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cle);
title('W_{NUMERIC}-W_{CURV}');

colormap(cptcmap('cool-warm'))
set(gcf, 'Color', 'w')
%%
tim = 400;
subplot(3,2,1)
% pcolor(X./1000, Y./1000, squeeze(sum(U(:,:,:, tim)-repmat(U(:,:,end,tim), [1 1 nz 1]), 3).*dz).');
pcolor(X./1000, Y./1000, squeeze(sum(Up(:,:,:, tim), 3).*dz).');

shading interp
colorbar;

subplot(3,2,2)
pcolor(X./1000, Y./1000, squeeze(sum(Vp(:,:,:, tim), 3).*dz).');
shading interp
colorbar;
% set(gca, 'clim', [-1.1 -0.9]);

subplot(3,2,3)
pcolor(X./1000, Y./1000, squeeze(UTRANS(:,:,1,tim)).');
shading interp
colorbar;

subplot(3,2,4)
pcolor(X./1000, Y./1000, squeeze(VTRANS(:,:,1,tim)).');
shading interp
colorbar;

%%
cl = [-1 1].*epsilon;

subplot(1,3,1)
% pcolor(X./1000, Y./1000, squeeze(sum(U(:,:,:, tim)-repmat(U(:,:,end,tim), [1 1 nz 1]), 3).*dz).');
contourf(X./1000, Y./1000, squeeze(ZETA(:,:)).'./(f));
set(gca, 'clim', cl);
shading interp
colorbar;

subplot(1,3,2)
% pcolor(X./1000, Y./1000, squeeze(sum(U(:,:,:, tim)-repmat(U(:,:,end,tim), [1 1 nz 1]), 3).*dz).');
contourf(X./1000, Y./1000, squeeze(ZETA(:,:)-omega(:,:)).'./(f));
set(gca, 'clim', cl);

shading interp
colorbar;

subplot(1,3,3)
% pcolor(X./1000, Y./1000, squeeze(sum(U(:,:,:, tim)-repmat(U(:,:,end,tim), [1 1 nz 1]), 3).*dz).');
contourf(X./1000, Y./1000, squeeze(omega(:,:)).'./(f));
set(gca, 'clim', cl);

shading interp
colorbar;
%%
zetas = ZETA - omega;

[zetasy, zetasx] = gradient(zetas, dy);
[omegay, omegax] = gradient(omega, dy);

we = (zetasx + omegay).*L/(f);
figure
pcolor(X./1000, Y./1000, squeeze(we).');
shading interp
colorbar

%%
r = Y./1000 - Y(floor(nx/2)-10)./1000;
zetas = ZETA-2*omega;
[var, ~] = gradient(ZETA, dy);
figure
plot(r, abs(epsilon.*Tavg(1,1)./(rho*f).*var(floor(nx/2),:)./f)*sqrt(2));%, r, abs(we(floor(nx/2),:)));
hold on
we = abs(WDIA) - abs(WSTERN);
plot(r,abs(we(floor(nx/2),:)));
hold off
grid on