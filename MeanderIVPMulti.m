%%% MeanderIVPMulti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
xfact = 120000;
deltax = 2000./xfact;
x = (0:deltax:300).*xfact;
f = 1e-4;
taumag =0.1/1030;
epsilon = .5;
ubarc = .5;
L = ubarc./(f*epsilon);
deltay = 500;
ys = -2*L:deltay:2*L;
ubart = ubarc.*exp(-(ys./L).^2/2);
A = 15e3;

Mutot = NaN(length(ys), length(x));
Mvtot = Mutot;
ytot = Mutot;
for ind=1:length(ys)
ubar = ubart(ind);
dudn = -ubar.*ys(ind)./(L.^2); 
epsk = A.*epsilon.*L.*(2*pi/xfact).^2;

rampwidth = 10*xfact;
xc = 275.*xfact;
facAmp = 1/2.*(1+tanh( (x-xc)./rampwidth));

% facAmp = facAmp.*exp(-( ( x-xc)./(6*xfact)).^2);
x1=x;
offset = ys(ind);
y = facAmp.*A.*sin(2*pi*x./(xfact))+offset;
ytot(ind,:) = y;

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

rampwidth = 40*xfact;
xc = 60.*xfact;
facTau = 1/2.*(1+tanh( (x-xc)./rampwidth));
% facTau = 1;
tau = facTau.*taumag.*ones(size(x));
taus = dot([real(tau); imag(tau)], [real(frntvec); imag(frntvec)]);
taun = dot([real(tau); imag(tau)], [-imag(frntvec); real(frntvec)]);

%%
% SOLVE ODES
omega = ubar*k;
zeta = -dudn + omega;

l = abs(cumtrapz(x1, abs(vels))); % int_x sqrt( 1+(dy/dx)^2)
% l = cumtrapz(real(velfr));
guess = [0 0]; %IVP bc u, v
out = meanderFrontODEIVP(l, omega, zeta, ubar, taus, taun, f, guess);
Mu = real(out.u.*frntvec + out.v.*1i.*frntvec);
Mv =  imag(out.u.*frntvec + out.v.*1i.*frntvec);

Mutot(ind,:) = Mu;
Mvtot(ind,:) = Mv;
end
%%
npaths = length(ys); nx = length(x);
xi = repmat(x, [1, length(ys)]).';
yi = reshape(ytot.', npaths*nx,1);
[X, Y] = meshgrid(x,ys);
muvec = reshape(Mutot.', npaths*nx,1);
mvvec = reshape(Mvtot.', 1,npaths*nx);
MUBAR = griddata(xi, yi, muvec, X, Y);
MVBAR = griddata(xi, yi, mvvec, X, Y);
[MUx, MUy] = gradient(MUBAR, deltax.*xfact);
[MVx, MVy] = gradient(MVBAR, deltay);
W = MUx+MVy;

ubfull = repmat(ubart, [nx, 1]);
ubfull = reshape(ubfull, 1, npaths*nx);
UBAR = griddata(xi, yi, ubfull, X, Y);

%%
inds = 15750:17500;
inds = 17200:17500;
uvecs = [.8 .9 .99].*ubarc;

subplot(3,1,1)
pcolor(x(inds), ys, MUBAR(:,inds)); shading interp
hold on
contour(x(inds), ys, UBAR(:,inds), uvecs, 'LineColor', 'k');
hold off
set(gca, 'clim', [-1 1]*1e-4*2000);
% set(gca, 'ylim', [-1 1]*30e3);
colormap(flipud(othercolor('RdBu11')));
colorbar;
title('U');

subplot(3,1,2)
pcolor(x(inds), ys, MVBAR(:,inds)); shading interp
hold on
contour(x(inds), ys, UBAR(:,inds), uvecs, 'LineColor', 'k');
hold off
set(gca, 'clim', [-1 1]*1e-4*20000);
% set(gca, 'ylim', [-1 1]*30e3);
colormap(flipud(othercolor('RdBu11')));
colorbar
title('V');

subplot(3,1,3)
wconts = linspace(-2e-4, 2e-4, 5);
pcolor(x(inds), ys, W(:,inds)); shading interp
% [c, h] = contourf(x(inds), ys, W(:,inds), wconts); 
% set(h, 'edgecolor', 'none');
hold on
contour(x(inds), ys, UBAR(:,inds), uvecs, 'LineColor', 'k');
hold off
set(gca, 'clim', [-1 1]*1e-4*2);
% set(gca, 'ylim', [-1 1]*30e3);
colormap(flipud(othercolor('RdBu11')));
colorbar
title('W');