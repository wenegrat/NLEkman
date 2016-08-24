function out = calcCircW(epsilon)

f = 1e-4;
tau = .1/1035;
thetas = 0:.01:2*pi;
cr = 60e3;
r = (0.01:.01:1).*2*cr;

ubarmax = epsilon.*f.*cr;
jetwidth = cr.*.6;

velstruct = r.*2./(jetwidth^2).*exp( - (r/jetwidth).^2);
velstruct = velstruct/max(abs(velstruct));
velstruct = ubarmax.*velstruct;
dudr = (1-2*r.^2./jetwidth.^2).*velstruct./r;

x = NaN(length(r), length(thetas));
y = x;
Myc = x;
Mx = x;
My = x;
umag = x;
maxom = 0;
for i=1:length(r)
   omega = velstruct(i)./r(i);
   if (abs(omega)>maxom); maxom = abs(omega); end
   zeta =  dudr(i) + omega;
   
   denom = (f+2*omega).*(f+zeta)-omega.^2;
%    Mr(i,:) = - (f+3*omega)./denom.*sin(thetas);
%    Mthet
   Mx(i,:) = (zeta - 2*omega)./denom.*sin(thetas).*cos(thetas).*tau;
   My(i,:) = -(f + omega + 2.*omega.*sin(thetas).^2 + zeta.*cos(thetas).^2)./denom.*tau;
   
   Myc(i,:) = -tau./(f+zeta).*ones(size(thetas));
   
   [x(i,:) y(i,:)] = pol2cart(thetas, r(i));
   uc(i,:) = velstruct(i).*(-sin(thetas));
   vc(i,:) = velstruct(i).*cos(thetas);
   umag(i,:) = velstruct(i).*ones(size(thetas));
end

nr = length(r); nth = length(thetas);

xvec = reshape(x, nr*nth, 1); 
yvec = reshape(y, nr*nth, 1);
Mxvec = reshape(Mx, nr*nth, 1);
Myvec = reshape(My, nr*nth, 1);

x = (-1:.01:1).*r(end); y = x;
deltax = x(2)-x(1); deltay=deltax;
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

umagvec = reshape(umag, nr*nth, 1);
UMAG = griddata(xvec, yvec, umagvec, X, Y);

out.x = X; out.y = Y;
out.UMAG = UMAG;
out.MU = MUBAR;
out.MV = MVBAR;
out.W = W;
out.WClassic = Wclassic;
end