deltax = 250;
x=0:deltax:1000e3;
f = 1e-4;
ubar = .25;
xc = floor(length(x)/2);
xc = x(xc);
amp = 12e3;
% amp = 6e3;

width = amp/1.5;
y = amp*exp( - ((x-xc)/width).^2);


% Try Cauchy distribution
% gamma = width;
% y =  (1./(pi*gamma* ( 1+ ( (x-xc)./gamma).^2)));
% y= y./max(y);
% y = amp.*y;

plot(x, y);
xc = 120e3;
rampwidth=80e3;
facTau = 1/2.*(1+tanh( (x-xc)./rampwidth));
% facTau = 1;
% tau = facTau.*taumag.*ones(size(x));
tau = .1.*facTau.*ones(size(x))/1035;
tau = tau-tau(1);


epsilon = .5;
ubarc = .5;
L = ubarc./(f*epsilon);
deltay = 500;
offsets = -amp/2:deltay:2*amp;
ubart = ubarc.*exp(-(offsets./L).^2/2);

%Allocate
u = NaN(length(offsets), length(x));
v = u; l = u; zeta = u; k = u; taus = u; taun = u; positions=u;ux=u; vx=u; 
svec = u; ubarall = u;

for i=1:length(offsets)
    disp([num2str(i),'/', num2str(length(offsets))]);
    off=offsets(i);
    pos = makeOffsetCurve(y, off, x);
    if isfinite(pos); 
    ubar = ubart(i);
    dudn = -ubar.*offsets(i)./(L.^2); 
    out = SolveFrontEkman(x, pos, ubar, dudn, tau, f);
    
    u(i,:) = out.u;
    ux(i,:) = out.ux;
    v(i,:) = out.v;
    vx(i,:) = out.vy;
    l(i,:) = out.l;
    zeta(i,:) = out.zeta;
    k(i,:) = out.omega./ubar;
    taus(i,:) = out.taus;
    taun(i,:) = out.taun;
    positions(i,:) = pos;
    svec(i,:) = out.tangent;
    ubarall(i,:) = ubar.*ones(size(x));
    end
end

%%
xi = repmat(x, [1, length(offsets)]).';
noff = length(offsets); nx = length(x);
lvec = reshape(l, noff*nx,1);
yi = reshape(positions.', noff*nx,1);
uvec = reshape(ux.', noff*nx, 1);
vvec = reshape(vx.', noff*nx, 1);
uevec = reshape(u.', noff.*nx, 1);
vevec = reshape(v.', noff*nx,1);
svecvec = reshape(svec.', noff*nx, 1);
kvec = reshape(k.', noff*nx,1);
ubarvec = reshape(ubarall.', noff*nx, 1);
zetavec = reshape(zeta.', noff*nx, 1);

[X, Y] = meshgrid(linspace(x(1), x(end), 4000), offsets);
deltax = X(1,2) - X(1,1);
deltay = Y(2,1) - Y(1,1); %% CHECK THIS!!!!

masky = isfinite(yi);
U = griddata(xi(masky), yi(masky), uvec(masky), X, Y);
V = griddata(xi(masky), yi(masky), vvec(masky), X, Y);
SVECX = griddata(xi(masky),yi(masky), svecvec(masky), X, Y);
K = griddata(xi(masky), yi(masky), kvec(masky), X, Y);
UBAR = griddata(xi(masky), yi(masky), ubarvec(masky), X, Y);
ZETABAR = griddata(xi(masky), yi(masky), zetavec(masky), X, Y);
Ue = griddata(xi(masky), yi(masky), uevec(masky), X, Y);
Ve = griddata(xi(masky), yi(masky), vevec(masky), X, Y);

mask = ones(size(Y)); 
ytemp = interp1(x, y, X(1,:));
for i=1:2000
    mask(Y(:,i)> ytemp(i)+amp, i) = NaN;
    mask(Y(:,i)< ytemp(i)+offsets(1), i) = NaN;
end
[Ux, Uy] = gradient(U, deltax, deltay);
[Vx, Vy] = gradient(V, deltax, deltay);

Ux = Ux.*mask; Uy = Uy.*mask; Vx=Vx.*mask; Vy=Vy.*mask;
WE = Ux + Vy; ZETABAR  = ZETABAR.*mask;

VortCart = Vx-Uy;


dvds = real(SVECX).*Vx + imag(SVECX).*Vy;
dudn = -imag(SVECX).*Ux + real(SVECX).*Uy;
uk = (real(SVECX).*U + imag(SVECX).*V).*K;
VortNat = dvds - dudn + uk;

%% CALCULATE TERMS IN BUDGET

%Along front advection of Ekman vorticity
[ZetaX, ZetaY] = gradient(VortCart, deltax, deltay);
ADV = UBAR.*( real(SVECX).*ZetaX + imag(SVECX).*ZetaY);
% ADV = ( ZetaY);
STRETCH = -(f+ZETABAR).*WE; % XX-REVISIT SIGN CONVENTION FOR UPWELLING

[ZBarX, ZBarY] = gradient(ZETABAR, deltax, deltay);
GRAD = -Ue.*( real(SVECX).*ZBarX + imag(SVECX).*ZBarY) - Ve.*(-imag(SVECX).*ZBarX + real(SVECX).*ZBarY);


%% CONFIRM VORTICITY
cl = [-1 1].*1e-7;
figure
subplot(2,1,1)
pcolor(X,Y, ADV); shading interp; colorbar;
set(gca, 'clim', cl);
    set(gca, 'xlim', [4.4e5 5.6e5]);
hold on
    plot(x, positions(ind,:))
hold off
subplot(2,1,2)
pcolor(X,Y, VortNat); shading interp; colorbar;
set(gca, 'clim', cl);
    set(gca, 'xlim', [4.4e5 5.6e5]);


%% CONTOURS
cl = [-1 1].*5e-8;

mask = ones(size(Y)); 
for i=1:length(x)
%     mask(Y(:,i)> y(i)+offsets(end)/4, i) = NaN;
%     mask(Y(:,i)< y(i)+offsets(1)/1, i) = NaN;
end
conts = linspace(cl(1), cl(end), 10);
figure
for i=1:5
    switch i
        case 1
            VAR = VortCart.*1e-4;
        case 2
            VAR = ADV;
        case 3
            VAR = STRETCH;
        case 4
            VAR = GRAD;
        case 5
            VAR = ADV-STRETCH-GRAD;
    end
   subplot(5,1,i) 
    [c, h] = contourf(X, Y, VAR.*mask, conts); shading interp
    set(h, 'edgecolor','none')
    colorbar;
    set(gca, 'clim', cl);

    set(gca, 'xlim', [4.8e5 5.8e5]);
%     set(gca, 'ylim', [offsets(1) offsets(end)]);
        hold on
    ind = 11;
    plot(x, positions(ind,:), 'k'); 
%     quiver(x, positions(ind,:), ux(ind,:), vx(ind,:), .25);
    hold off
%     axis equal
end

%%
cl = [-1 1].*1e-7;

mask = ones(size(Y)); 
for i=1:length(x)
    mask(Y(:,i)> y(i)+offsets(end)/2, i) = NaN;
    mask(Y(:,i)< y(i)+offsets(1)/1, i) = NaN;
end
conts = linspace(cl(1), cl(end), 10);

    [c, h] = contourf(X, Y, VortCart.*1e-4, conts); shading interp
    set(h, 'edgecolor','none')
    colorbar;
    set(gca, 'clim', cl);
        hold on
    ind = 9;
    plot(x, positions(ind,:), 'k'); 
    quiver(x, positions(ind,:), ux(ind,:), vx(ind,:));
    hold off
    set(gca, 'xlim', [4.4e5 6e5]);
%     set(gca, 'ylim', [offsets(1) offsets(end)]);

    axis equal


%%
% LINE PLOT VERSION
xnorm = 2*pi*ubarc./(f);
lxpos = .125;
fs = 12;
xl = [4.24e5 8e5];
index = find(x>xl(1), 1)-1;
indexe = find(x>xl(2), 1);
xl = [0 6];
xn = x./xnorm;
xn = xn - xn(index);
gap = [.02 .05]; margh = .1; margw=.1;
figure
subtightplot(3,1,1,gap, margh, margw)
ind = 13; %11 works well.
ADVI = interp2(X, Y, ADV, x, positions(ind,:));
STRETCHI = interp2(X, Y, STRETCH, x, positions(ind,:));
GRADI = interp2(X,Y, GRAD, x, positions(ind,:));

% PLOTS
rp = 6;
plot(xn, positions(ind,:)./xnorm, 'LineWidth', 2);
hold on
set(gca, 'ColorOrderIndex', 1);
q = quiver(xn(1:rp:end), positions(ind,1:rp:end)./xnorm, ux(ind,1:rp:end), vx(ind,1:rp:end));
set(q, 'ShowArrowHead', 'off');
hold off
axis equal
set(gca, 'xlim', xl);
grid on
% title(num2str(Y(ind,1)));
% axis equal
set(gca, 'ylim', [-1.5 .5]);
set(gca, 'XTickLabel', []);
ylabel('$\hat{y}$', 'Interpreter', 'Latex');
t = text(lxpos, 0.275, 'Ekman Transport');
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
vortnorm = ubarc.*f.*20/L;
set(gca, 'FontSize', fs);

% vortnorm = taumag./(f.*L).*ubarc;
subtightplot(3,1,2,gap, margh, margw)
plot(xn, ADVI./vortnorm, 'LineWidth', 2);
hold on
plot(xn, STRETCHI./vortnorm, 'LineWidth', 2);
plot(xn, GRADI./vortnorm, 'LineWidth', 2);
hold off
set(gca, 'xlim', xl);
grid on
% set(gca,'ylim', [-10 20]);
% set(gca, 'ylim', [-1 1]);
legend('ADV', 'STRETCH', 'GRAD', 'Location', 'NorthEast');
set(gca, 'XTickLabel', []);
t = text(lxpos, .80, 'Vorticity Budget');
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
set(gca, 'FontSize', fs);

% ylabel('\hat
vortnorm = vortnorm.*L;
subtightplot(3,1,3,gap, margh, margw)
plot(xn(index:indexe), cumtrapz(x(index:indexe), ADVI(index:indexe))./vortnorm, 'LineWidth', 2);
hold on
plot(xn(index:indexe), cumtrapz(x(index:indexe), STRETCHI(index:indexe))./vortnorm, 'LineWidth', 2);
plot(xn(index:indexe), cumtrapz(x(index:indexe), GRADI(index:indexe)./vortnorm), 'LineWidth', 2);
% plot(xn(index:indexe), cumtrapz(x(index:indexe), ADVI(index:indexe) - ...
%     STRETCHI(index:indexe)-GRADI(index:indexe))./vortnorm, 'LineWidth', 2, 'LineStyle', '--');

hold off
set(gca, 'xlim', xl);
grid on
xlabel('$\hat{x}$', 'Interpreter', 'Latex');
t = text(lxpos, .75, 'Cumulative Budget');
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
set(gca, 'FontSize', fs);
set(gcf,'Color', 'w', 'Position', [665   376   548   592])

% TO EXPORT, FIX LEGEND POSITION THEN
% export_fig('VorticityBudgetNew.eps', '-eps', '-painters', '-q101', '-p01')
%%
figure
subplot(2,1,1);
xc = floor(length(x)/2);
xc = x(xc);
inds =440:575;
xn = x(inds)./nrm; yn = y(inds)/nrm; 
quiver(xn-xn(1), yn, out.ux(inds), out.vy(inds));
hold on; plot(xn-xn(1), yn, 'LineWidth', 2); hold off
width = 2*pi./(f./ubar);
axis equal
set(gca, 'ylim', [-1 1].*1);
set(gca, 'xlim', [0 7.5])
hold on
% quiver(.5, .75, tau(end), 0)
hold off
grid on
xlabel('$x/\lambda_n$', 'Interpreter', 'Latex');
ylabel('$y/\lambda_n$', 'Interpreter', 'Latex');
set(gca, 'FontSize', 20);

subplot(2,1,2)
% plot(xn-xn(1), vort(inds), 'LineWidth', 2);
hold on
plot(xn-xn(1), adv(inds), 'LineWidth', 2);
plot(xn-xn(1), tilting(inds), 'LineWidth', 2);
plot(xn-xn(1), ekadv(inds), 'LineWidth', 2);
plot(xn-xn(1), forcing(inds), 'LineWidth', 2);
legend('Advection','Stretching','Parametric', 'Forcing', 'location', 'NorthWest');
hold off
set(gca, 'xlim', [0 7.5])
grid on
box on
set(gca, 'FontSize', 20);
xlabel('$x/\lambda_n$', 'Interpreter', 'Latex');
ylabel('\int \zeta dz (m/s)');
set(gcf, 'Color', 'w', 'Position', [   675   374   805   601]);

%%
fs=14;
figure
subplot(3,1,1);
xc = floor(length(x)/2);
xc = x(xc);
inds =440:575;
xn = x(inds)./nrm; yn = y(inds)/nrm; 
quiver(xn-xn(1), yn, out.ux(inds), out.vy(inds));
hold on; plot(xn-xn(1), yn, 'LineWidth', 2); hold off
width = 2*pi./(f./ubar);
axis equal
set(gca, 'ylim', [-1 1].*1);
set(gca, 'xlim', [0 7.5])
% set(gca, 'xlim', [xc-4*width xc+4*width]);
% set(gca, 'ylim', [-2.5*amp 2.5*amp])
% cb = colorbar;
% set(cb, 'visible', 'off')
% axis equal
hold on
% quiver(.5, .75, tau(end), 0)
hold off
grid on
xlabel('$x/\lambda_n$', 'Interpreter', 'Latex');
ylabel('$y/\lambda_n$', 'Interpreter', 'Latex');
set(gca, 'FontSize', fs);
title('Ekman Transport');

subplot(3,1,2)
% plot(xn-xn(1), vort(inds), 'LineWidth', 2);
hold on
plot(xn-xn(1), adv(inds), 'LineWidth', 2);
plot(xn-xn(1), tilting(inds), 'LineWidth', 2);
plot(xn-xn(1), ekadv(inds), 'LineWidth', 2);
plot(xn-xn(1), forcing(inds), 'LineWidth', 2);
legend('Advection','Stretching','Parametric', 'Forcing', 'location', 'NorthWest');
hold off
set(gca, 'xlim', [0 7.5])
grid on
box on
set(gca, 'FontSize', fs);
xlabel('$x/\lambda_n$', 'Interpreter', 'Latex');
title('Vorticity Budget Terms');
ylabel('m/s^2');

subplot(3,1,3)
plot(xn-xn(1), cumtrapz(out.l(inds),adv(inds)./ubar), 'LineWidth', 2);
hold on
plot(xn-xn(1), cumtrapz(out.l(inds),tilting(inds)./ubar), 'LineWidth', 2);
plot(xn-xn(1), cumtrapz(out.l(inds),ekadv(inds)./ubar), 'LineWidth', 2);
plot(xn-xn(1), cumtrapz(out.l(inds),forcing(inds)./ubar), 'LineWidth', 2);
hold off
set(gca, 'xlim', [0 7.5])
grid on
box on
set(gca, 'FontSize', fs);
xlabel('$x/\lambda_n$', 'Interpreter', 'Latex');
title('Cumulative Vorticity Budget');
ylabel('m^2/s^2')
set(gcf, 'Color', 'w', 'Position', [  269   287   623   685]);