%% COMPARE CARTESIAN AND CURVATURE VORTICITY COMPONENTS
thetas = 0:.01:2*pi;
cr = 60e3;
r = (0.01:.01:1).*2*cr;

ubarmax = 1;
jetwidth = cr.*.6;

%Gaussian Vel
velstruct = ubarmax.*exp(- ((r-cr)/jetwidth).^2);
dudr = -(r-cr).*2./(jetwidth).^2.*velstruct;

fullVel = NaN(length(r), length(thetas));
zetas = fullVel;
zetac = fullVel;
x = fullVel;
y = fullVel;
for i =1:length(r)
    ubar = velstruct(i);
    fullVel(i,:) = -ubar.*sin(thetas) + 1i*ubar.*cos(thetas);
    zetas(i,:) = dudr(i);
    zetac(i,:) = ubar./r(i);
    
    [x(i,:), y(i,:)] = pol2cart(thetas, r(i));
end

nr = length(r); nth = length(thetas);

xvec = reshape(x, nr*nth, 1); 
yvec = reshape(y, nr*nth, 1);
zetasvec = reshape(zetas, nr*nth, 1);
zetacvec = reshape(zetac, nr*nth, 1);
velvec = reshape(fullVel, nr*nth, 1);

x = (-1:.01:1).*r(end); y = (-1:.005:1).*r(end);
deltax = x(2)-x(1); deltay=y(2)-y(1);
[X, Y] = meshgrid(x, y);
U = griddata(xvec, yvec, real(velvec), X, Y);
V = griddata(xvec, yvec, imag(velvec), X, Y);
[Ux, Uy] = gradient(U, deltax, deltay);
[Vx, Vy] = gradient(V, deltax, deltay);
ZETACART = Vx - Uy;
DIV = Ux + Vy;

ZETASH = griddata(xvec, yvec, zetasvec, X, Y);
ZETACU = griddata(xvec, yvec, zetacvec, X, Y);

figure
subplot(2,2,1)
pcolor(x, y, Vx); shading interp
set(gca, 'clim', [-2 10].*1e-5)
subplot(2,2,2)
pcolor(x, y, -Uy); shading interp
set(gca, 'clim', [-2 10].*1e-5)
subplot(2,2,3)
pcolor(x, y, ZETASH); shading interp
set(gca, 'clim', [-2 10].*1e-5)
subplot(2,2,4)
pcolor(x, y, ZETACU); shading interp
set(gca, 'clim', [-2 10].*1e-5)
%%
figure
subplot(2,1,1)
pcolor(ZETACART) ;shading interp 
set(gca, 'clim', [-2 10].*1e-5)
colorbar;
subplot(2,1,2)
pcolor(ZETASH+ZETACU); shading interp
colorbar;
set(gca, 'clim', [-2 10].*1e-5)
