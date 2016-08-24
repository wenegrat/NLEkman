% SOLVE VIA IVP

xfact = 120000;
deltax = 2000./xfact;
x = (0:deltax:600).*xfact;
f = 1e-4;
taumag =0.1/1030;
epsilon = .5;
ubarc =1;
L = ubarc./(f*epsilon);
deltay = 500;
ys = -2*L:deltay:2*L;
ubart = ubarc.*exp(-(ys./L).^2/2);
A = 12e3;

ind = 49; % find from MeanderFront
ubar = ubart(ind);
ubar = ubarc;
dudn = -ubar.*ys(ind)./(L.^2); 
epsk = A.*epsilon.*L.*(2*pi/xfact).^2;

rampwidth = 10*xfact;
xc = 200.*xfact;
facAmp = 1/2.*(1+tanh( (x-xc)./rampwidth));
facAmp=1;
% facAmp = facAmp.*sin(2*pi*x./(12*xfact));
% facAmp = facAmp.*exp(-( ( x-xc)./(6*xfact)).^2);
x1=x;
y = facAmp.*A.*sin(2*pi*x./(xfact));

% Define alternate style of curvature.
% vel = ubar;
% maxa = 180*pi./180;
% angle = maxa.*cos((2*pi./(3*xfact)).*x);
% velfr = vel.*exp(1i*angle);
% y = cumtrapz(x, imag(velfr));
% x = cumtrapz(x,real(velfr));
% y = facAmp.*A.*y./max(y);
% y = 2*A*sinc( (x-200*xfact)./(1/2.*xfact));
% x = facAmp.*A.*cos(2*pi*x./(xfact));

%Determine curvature
dx  = gradient(x, deltax*xfact);
ddx = gradient(dx, deltax*xfact);
dy  = gradient(y, deltax*xfact);
dy(1:end-1) = (y(2:end)-y(1:end-1))./(deltax.*xfact);
dy(end) = dy(1);
ddy = gradient(dy, deltax*xfact);
ddy(1:end-1) = (dy(2:end)-dy(1:end-1))./(deltax.*xfact);
ddy(end) =ddy(1);
num   = dx .* ddy - ddx .* dy;
denom = dx .* dx + dy .* dy;
denom = sqrt(denom);
denom = denom.* denom.* denom;
k = num ./ denom;
k(denom < 0) = NaN;
disp(['Eps-Omega = ', num2str(max(k).*ubar./f,2)]);

% Determine s and n vectors
du = gradient(x, deltax*xfact); %Temporary velocity gradients
dv = gradient(y, deltax*xfact); %
dv(1:end-1) = (y(2:end)-y(1:end-1))./(deltax*xfact);
dv(end)=dv(1);
vels = du+1i*dv;  %Tangent Vectors at each spot.
frntvec = vels./abs(vels); %Normalized;

rampwidth = 100*xfact;
xc = 200.*xfact;
facTau = 1/2.*(1+tanh( (x-xc)./rampwidth));
% facTau = 1;
tau = facTau.*taumag.*ones(size(x));
taus = dot([real(tau); imag(tau)], [real(frntvec); imag(frntvec)]);
taun = dot([real(tau); imag(tau)], [-imag(frntvec); real(frntvec)]);
taus = tau; taun=0*taun;
figure
subplot(3,1,1)
plot(x, y);
ylabel('Front position');
subplot(3,1,2)
plot(x, tau);
ylabel('\tau');
%%
% SOLVE ODES
omega = ubar*k;
zeta = -dudn + omega;

l = abs(cumtrapz(x1, abs(vels))); % int_x sqrt( 1+(dy/dx)^2)
% l = cumtrapz(real(velfr));
guess = [0 0]; %IVP bc u, v
% taus = mean(taus)+0*taus; taun = 0*taun;
out = meanderFrontODEIVP(l, omega, zeta, ubar, taus, taun, f, guess);
Mu = real(out.u.*frntvec + out.v.*1i.*frntvec);
Mv =  imag(out.u.*frntvec + out.v.*1i.*frntvec);

utheory = (f+zeta + omega)./((f+2*omega).*(f+zeta)).*taun;
vtheory = -(f+3*omega)./((f+2*omega).*(f+zeta)).*taus;
Mut = real(utheory.*frntvec + vtheory.*1i.*frntvec);
Mvt = imag(utheory.*frntvec + vtheory.*1i.*frntvec);
[b, a] = butter(4, 4*.025);
Mvs = filtfilt(b,a,Mv);
subplot(3,1,3)
plot(x,Mv);
hold on
plot(x, Mu);
hold off
ylabel('V');

%%
figure
subplot(3,1,1)
inds = length(x)-250:length(x)-3;
inds = length(x)-5000:length(x)-250;
inds = 11000:12500;
% inds = 1/deltax*195:1*205/deltax;
quiver(x(inds), y(inds), Mu(inds), Mv(inds));
hold on
plot(x(inds), y(inds), 'LineWidth',2);
hold off
axis equal
subplot(3,1,2)
quiver(x(inds), y(inds), Mut(inds), Mvt(inds));
hold on
plot(x(inds), y(inds), 'LineWidth',2);
hold off
axis equal
subplot(3,1,3)
plot(x(inds), Mv(inds));
hold on
plot(x(inds), Mvt(inds), 'LineWidth', 2);
plot(x(inds), Mvs(inds), 'LineWidth', 2);
% plot(x(inds), Mu(inds));
hold off
%%
F = sqrt((f+zeta).*(f+2*omega));
% k = F./ubar;
% 
% ei = 100;
% l = cumtrapz(

figure
subplot(4,1,1)
plotyy(l(inds), -(f+3*omega(inds))./(F(inds).^2), l(inds), (f+zeta(inds))./(F(inds).^2))
title('V factor     U factor');
subplot(4,1,2)
plot(l(inds), taus(inds));
hold on
plot(l(inds), taun(inds));
hold off
title('\tau_s     \tau_n');
subplot(4,1,3)
plot(l(inds), -(f+3*omega(inds))./(F(inds).^2).*taus(inds))
hold on;
plot( l(inds), out.v(inds))
title('V_{theory}    V_{numeric}');
subplot(4,1,4)
plotyy(l(inds), (f+zeta(inds))./(F(inds).^2).*taun(inds), l(inds), out.u(inds))
title('U_{theory}    U_{numeric}');
%%
trapz(l, F./ubar)./(2*pi*600)

inds = 11000:13500;
figure
subplot(3,1,1)
plot(l(inds), y(inds))
grid on
ylabel('Front Y Position');
% set(gca, 
subplot(3,1,2)
vfact = -(f+3*omega(inds))./(F(inds).^2);
plot(l(inds), vfact./mean(vfact));
hold on
ufact = (f+zeta(inds))./(F(inds).^2);
plot( l(inds), abs(ufact./mean(ufact(1))))
hold off
legend('1/F_V', '1/F_U');
grid on
ylabel('Normalized');

subplot(3,1,3)
plot(l(inds), taus(inds));
hold on
plot(l(inds), taun(inds));
hold off
legend('\tau_s', '\tau_n');
grid on
ylabel('\tau/\rho');