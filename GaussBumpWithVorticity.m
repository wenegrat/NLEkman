x=0:1e3:1000e3;
f = 1e-4;
ubar = .25;
xc = floor(length(x)/2);
xc = x(xc);
amp = 12e3;
width = amp/2;
y = amp*exp( - ((x-xc)/width).^2);

plot(x, y);
xc = 120e3;
rampwidth=80e3;
facTau = 1/2.*(1+tanh( (x-xc)./rampwidth));
% facTau = 1;
% tau = facTau.*taumag.*ones(size(x));

tau = .1.*facTau.*ones(size(x))/1035;
tau = tau-tau(1);
zetabase = 0;.1*f;
out = SolveFrontEkman(x, y, ubar, zetabase, tau, f);

% Vorticity Budget
vort = gradient(out.v.', out.l.');
% vort = (out.v(2:end)-out.v(1:end-1))./(out.l(2:end)-out.l(1:end-1));
adv = ubar.*4*del2(out.v, out.l);
tilting = (f+2*out.omega).*gradient(out.u, out.l);
ekadv = out.u.*gradient(2*out.omega, out.l);
forcing = gradient(out.taun, out.l);

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