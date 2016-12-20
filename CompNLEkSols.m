%%
% Compare theoretical solutions
%%
statefile = 'state.nc'; diagfile = 'diag.nc'; etanfile = 'etan.nc';
f = 1e-4;
rho = 1035;


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


U = GetVar(statefile, diagfile, {'UVEL', '(1)'}, {0, 0, 0, 0});
V = GetVar(statefile, diagfile, {'VVEL', '(1)'}, {0, 0, 0, 0});
W = GetVar(statefile, diagfile, {'WVEL', '(1)'}, {0, 0, 0, 0});
%%
TAU = GetVar(statefile, etanfile, {'oceTAUX', '(1)'}, {0, 0, [1 1], 0});
ETA = GetVar(statefile, etanfile, {'ETAN', '(1)'}, {0, 0, [1 1], 0});
%%
[omega, phi, UB, VB] = calculateOmega(ETA,U, V, dx, dy);

% ZETAF = DPeriodic(V, dx, 'x') - DPeriodic(U, dy, 'y');
%% TEMP
% TAU = .1.*ones(size(ZETA));
Up = U - repmat(U(:,:,end,:), [1 1 nz 1]);
Vp = V - repmat(V(:,:,end,:), [1 1 nz 1]);

Upt = sum(Up, 3).*dz;
Vpt = sum(Vp, 3).*dz;
WDIA = (DPeriodic(Upt, dx, 'x') + DPeriodic(Vpt, dy, 'y'));

%%
ZETA = DPeriodic(V(:,:,end,:), dx, 'x') - DPeriodic(U(:,:,end,:), dy, 'y');
%% STERN SOLUTIONS
VTRANS_STERN = -repmat(TAU,[1 1 1, 1])./(rho*(f+ZETA));
% VTRANS_STERN = -TAU./(rho*(f+ZETA));

WSTERN = gradient(VTRANS_STERN, dy);


%% CURVED SOLUTIONS
theta = phi-pi/2;

denom = ( (f + 2.*omega).*(f+ZETA) - omega.^2);
VTRANS = - (f + omega + 2.*omega.*sin(theta).^2 + ZETA.*cos(theta).^2)./denom.*TAU./rho;
UTRANS =  (ZETA-2.*omega)./denom.*TAU./rho.*sin(theta).*cos(theta);

[Vy, Vx, ~, ~] = gradient(VTRANS, dx);
[Uy, Ux, ~, ~] = gradient(UTRANS, dx);

WCURV = Ux + Vy;

%%
cl = [-1 1].*1e-6;
cle = cl./20;

dep = 1;
tim = 480;
figure
subplot(3,2,1)
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

subplot(3,2,3)
pcolor(X./1000, Y./1000, squeeze(WSTERN(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cl);
title('W_{STERN}');

subplot(3,2,4)
pcolor(X./1000, Y./1000, squeeze(WDIA(:,:,dep, tim) - WSTERN(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cle);
title('W_{NUMERIC}-W_{STERN}');

subplot(3,2,5)
pcolor(X./1000, Y./1000, squeeze(WCURV(:,:,1, tim)).');
shading interp
colorbar
set(gca, 'clim', cl);
title('W_{CURV}');

subplot(3,2,6)
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