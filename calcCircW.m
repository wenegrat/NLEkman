function out = calcCircW(epsilon)

f = 1e-4;
tau = .1/1035;
thetas = 0:.01:2*pi;
cr = 60e3;
r = (0.01:.01:1).*2*cr;

ubarmax = epsilon.*f.*cr;
jetwidth = cr.*.6;

% Gaussian SSH
velstruct = r.*2./(jetwidth^2).*exp( - (r/jetwidth).^2);
velstruct = velstruct/max(abs(velstruct));
velstruct = ubarmax.*velstruct;
dudr = (1-2*r.^2./jetwidth.^2).*velstruct./r;

%Gaussian Vel
% velstruct = ubarmax.*exp(- ((r-cr)/jetwidth).^2);
% dudr = -(r-cr).*2./(jetwidth).^2.*velstruct;

x = NaN(length(r), length(thetas));
y = x;
Myc = x;
Mx = x;
My = x;
umag = x;
zetaf=x;
omegaf=x;
Myeps = x;
Mxeps = x;
Mvclass = x;
Mvlin = x;
maxom = 0;
for i=1:length(r)
   omega = velstruct(i)./r(i);
   if (abs(omega)>maxom); maxom = abs(omega); end
   zeta =  dudr(i) + omega;
   
   denom = (f+2*omega).*(f+zeta)-omega.^2;

   Mx(i,:) = (zeta - 2*omega)./denom.*sin(thetas).*cos(thetas).*tau;
   My(i,:) = -(f + omega + 2.*omega.*sin(thetas).^2 + zeta.*cos(thetas).^2)./denom.*tau;
   
   Mr(i,:) = - (f + 3*omega)./denom.*sin(thetas).*tau;
   Mth(i,:) = -(f + omega+zeta)./denom.*cos(thetas).*tau;

   Mx(i,:) = Mr(i,:).*cos(thetas) - Mth(i,:).*sin(thetas);
   My(i,:) = Mr(i,:).*sin(thetas) + Mth(i,:).*cos(thetas);
   
   Myc(i,:) = -tau./(f+zeta).*ones(size(thetas));%% This is wrong, need to consider angle!
   
   Myeps(i,:) = -(f - (zeta-omega).*sin(thetas).^2 - omega.*cos(thetas).^2).*tau./f.^2;
   Mxeps(i,:) = (zeta-2*omega).*sin(thetas).*cos(thetas).*tau;
   
   Mvclass(i,:) = - tau./(f+zeta);
   Mvlin(i,:) = -tau./f;
   
   [x(i,:), y(i,:)] = pol2cart(thetas, r(i));
   uc(i,:) = velstruct(i).*(-sin(thetas));
   vc(i,:) = velstruct(i).*cos(thetas);
   umag(i,:) = velstruct(i).*ones(size(thetas));
   omegaf(i,:) = omega;
   zetaf(i,:) = zeta;
end

nr = length(r); nth = length(thetas);

xvec = reshape(x, nr*nth, 1); 
yvec = reshape(y, nr*nth, 1);
Mxvec = reshape(Mx, nr*nth, 1);
Myvec = reshape(My, nr*nth, 1);

x = (-1:.01:1).*r(end); y = (-1:.005:1).*r(end);
deltax = x(2)-x(1); deltay=y(2)-y(1);
[X, Y] = meshgrid(x, y);

MUBAR = griddata(xvec, yvec, Mxvec, X, Y);
MVBAR = griddata(xvec, yvec, Myvec, X, Y);
[MUx, MUy] = gradient(MUBAR, deltax);
[MVx, MVy] = gradient(MVBAR, deltay);
W = MUx+MVy;

Mycvec = reshape(Myc, nr*nth, 1);
MVBARC = griddata(xvec, yvec, Mycvec, X, Y);
[MVx, MVy] = gradient(MVBARC, deltay);
Wclassic = MVy;

omegavec = reshape(omegaf, nr*nth,1);
zetavec = reshape(zetaf, nr*nth,1);
OMEGA = griddata(xvec, yvec, omegavec, X, Y);
ZETA = griddata(xvec, yvec, zetavec, X, Y);

umagvec = reshape(umag, nr*nth, 1);
UMAG = griddata(xvec, yvec, umagvec, X, Y);

myepsvec = reshape(Myeps, nr*nth, 1);
MVEPS = griddata(xvec, yvec, myepsvec, X, Y);
mxepsvec = reshape(Mxeps, nr*nth, 1);
MUEPS = griddata(xvec, yvec, mxepsvec, X, Y);
mvclassvec = reshape(Mvclass,nr*nth,1);
MVCLASS = griddata(xvec, yvec, mvclassvec, X, Y);
mvlinvec = reshape(Mvlin, nr*nth,1);
MVLIN = griddata(xvec, yvec, mvlinvec, X, Y);

ucvec = reshape(uc, nr*nth, 1);
vcvec = reshape(vc, nr*nth,1);
UC = griddata(xvec, yvec, ucvec, X, Y);
VC = griddata(xvec, yvec, vcvec, X, Y);
[UX, UY] = gradient(UC, deltay);
[VX, VY] = gradient(VC, deltax);
ZETACART = VX-UY;
MVCLASS = -tau./(f+ZETACART);
% Wclassic = UX+VY;

out.x = x; out.y = y;
out.UMAG = UMAG;
out.MU = MUBAR;
out.MV = MVBAR;
out.W = W;
out.WClassic = Wclassic;
out.OMEGA = OMEGA;
out.ZETA = ZETA;
out.MVEPS = MVEPS;
out.MUEPS = MUEPS;
out.MVCLASS= MVCLASS;
out.MVLIN = MVLIN;
end
