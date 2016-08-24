%%
% SOLVE THE CURVED FRONT EKMAN TRANSPORT
%%
% close all
% PARAMETERS
% l(x,y)- Parameterized along front coordinate
% deltax = (.01);
xfact = 120000;
deltax = 500./xfact;
x = (0:deltax:100).*xfact;
% x = A*cos(2*pi
% x = (0:deltax:200).*xfact;
% ubar - Along front speed (uniform in l), with cross front variation

% f - Coriolis
f = 1e-4;
% tau
taumag =0.1/1030;

% DERIVED QUANTITIES
% epsilon
epsilon = .5;
% L = 7e3;
% ubarc = (f*L*epsilon);
ubarc = .5;
L = ubarc./(f*epsilon);
% L = ubarc./(f.*epsilon);
deltay = 500;
ys = -3*L:deltay:3*L;
ubart = ubarc.*exp(-(ys./L).^2/2);
% dubar/dn - shear vorticity
% dudn = -ubarc./L;

% epsilonr/epsilon = A/L
maxA = .5*(xfact./(2*pi)).^2./(L*epsilon);
% max
Avecs = [0 4e3 8e3 12e3];
if Avecs(end)>maxA
    disp('Warning: A too large');
end
epsk = Avecs.*epsilon.*L.*(2*pi/xfact).^2;
er = epsk./epsilon;
% figure 
%Note by Holton's discussion, requires that:
% f*ubar_g < f^2./(4*k)
% What is a meaningful definition of ubar_g here? (probably ubarc (ie
% straight frontal case).
% So ubarc has to be small enough, and K small enough, k max =
% A*(2*pi./xfact)^2.
A_max = (xfact./(2*pi)).^2.*f./(4*ubarc);

for i=4;1:4;
% i =2;
A =Avecs(i);
ytot = NaN(length(ys), length(x));
Mutot = ytot;
Mvtot = ytot;
dudntot = NaN(length(ys),1);
zetatot = Mvtot;
omegatot = Mvtot;
ktot = Mvtot;
kf = Mvtot;
for j=  75;1:length(ubart);
ubar = ubart(j);
dudn = -ubar.*ys(j)./(L.^2); % IS THIS CORRECT FOR GAUSSIAN??
dudntot(j) = dudn;
offset = ys(j);
rampwidth = 1e6;
xc = x(end)/3;
fac = 1/2.*(1+tanh( (x-xc)./rampwidth));
% fac = sin(2*pi*x./(20*xfact));
% fac = ones(size(x));
y = fac.*A.*sin(2*pi*x./(xfact))+offset;
% y = A*exp(-((x-xfact/2)/(2*A)).^2)+offset;
% y(1:481) = offset;
ytot(j,:) = y;
x1 = x;
% x = A*cos(2*pi*x./xfact);
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
ktot(j,:) = k;
% k(num==0) = NaN;

du = gradient(x, deltax*xfact); %Temporary velocity gradients
dv = gradient(y, deltax*xfact); %
dv(1:end-1) = (y(2:end)-y(1:end-1))./(deltax*xfact);
dv(end)=dv(1);

vels = du+1i*dv;  %Tangent Vectors at each spot.
frntvec = vels./abs(vels); %Normalized;

tau = taumag.*ones(size(x));
taus = dot([real(tau); imag(tau)], [real(frntvec); imag(frntvec)]);
% taus = mean(taus) + zeros(size(taus));
taun = dot([real(tau); imag(tau)], [-imag(frntvec); real(frntvec)]);
% taun  = mean(taun) + zeros(size(taun));
%%%%
% SOLVE ODES
omega = ubar*k;
zeta = -dudn + omega;
zetatot(j,:) = zeta;
omegatot(j,:) = omega;
l = abs(cumtrapz(x, vels)); % int_x sqrt( 1+(dy/dx)^2)
l = abs(cumtrapz(x1, abs(vels))); % int_x sqrt( 1+(dy/dx)^2)
% l = cumtrapz(abs(ubar.*(frntvec))).*xfact;
if j==1
    guessvec = [taun./(f+2*omega); -taus./f];
else
    guessvec = [taun./(f+2*omega); -taus./(f+zeta)];
end
out = meanderFrontODEs(l, omega, zeta, ubar, taus, taun, f, guessvec);
Mu = real(out.u.*frntvec + out.v.*1i.*frntvec);
Mv =  imag(out.u.*frntvec + out.v.*1i.*frntvec);
% Mv = imag(out.v.*(1i*frntvec));

utheory = (f+zeta + omega)./((f+2*omega).*(f+zeta)).*taun;
vtheory = -(f+3*omega)./((f+2*omega).*(f+zeta)).*taus;
Mut = real(utheory.*frntvec + vtheory.*1i.*frntvec);
Mvt = imag(utheory.*frntvec + vtheory.*1i.*frntvec);

Mutot(j,:) = Mu;
Mvtot(j,:) = Mv;


 ktau = gradient(angle(frntvec),l);
 kf(j,:) = sqrt((f+zeta).*(f+2*omega))./ubar;
end


%%
% UBAR = pathToGrid(repmat(ubart.', [1 101]), ytot, x,ys);

% MUBAR = pathToGrid(Mutot, ytot, x, ys);

npaths = length(ys); nx = length(x);
xi = repmat(x, [1, length(ys)]).';
yi = reshape(ytot.', npaths*nx,1);
% xi = repmat(x, [1, length(ys)-2]).';
% ytottemp = ytot;
% ytottemp(65,:) = [];
% ytottemp(28,:) =[];
% yi = reshape(ytottemp.', (npaths-2)*nx,1);
% Mutottemp = Mutot;
% Mutottemp(65,:) = [];
% Mutottemp(28,:) = [];
% Mvtottemp = Mvtot;
% Mvtottemp(65,:) = [];
% Mvtottemp(28,:) = [];
% muvec = reshape(Mutottemp.', (npaths-2)*nx,1);
% mvvec = reshape(Mvtottemp.', 1,(npaths-2)*nx);

muvec = reshape(Mutot.', npaths*nx,1);
mvvec = reshape(Mvtot.', 1,npaths*nx);

ubfull = repmat(ubart, [nx, 1]);
ubfull = reshape(ubfull, 1, npaths*nx);
[X, Y] = meshgrid(x,ys);
MUBAR = griddata(xi, yi, muvec, X, Y);
MVBAR = griddata(xi, yi, mvvec, X, Y);
UBAR = griddata(xi, yi, ubfull, X, Y);
% MVBAR = pathToGrid(Mvtot, ytot, x, ys);
[MUx, MUy] = gradient(MUBAR, deltax.*xfact);
[MVx, MVy] = gradient(MVBAR, deltay);
W = MUx+MVy;

ZETABAR = griddata(xi, yi, reshape(zetatot.', 1, npaths*nx), X, Y);
OMEGABAR = griddata(xi,yi, reshape(omegatot.', 1, npaths*nx), X, Y);
%%
uvecs = [.8 .9 .99].*ubarc;
subplot(4,1,i)
pcolor(x, ys, W); shading interp
hold on
% contour(x, ys, UBAR,  uvecs,'LineColor', 'k');
hold off
set(gca, 'clim', [-1 1]*1e-4);
set(gca, 'ylim', [-1 1]*30e3);
colormap(flipud(othercolor('RdBu11')));
axis equal
title(['Amp: ', num2str(A/1000), '  \epsilon_K =', num2str(epsk(i),2), '  \epsilon = ', num2str(epsilon,2), '  u = ',num2str(ubarc,2), '  \epsilon_k/\epsilon = ', num2str(epsk(i)./epsilon, 2)]);
end