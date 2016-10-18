function out = calcCircW(epsilon, tau, f, cr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CalcCircW
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function calculates the NL Ekman velocities at a circular eddy
%
% Inputs: 
%       epsilon - ubar/(fL)
%       tau     - tau_o/rho (uniform tau)
%       f       - Coriolis
%       cr      - Radius in meters
%
%       out = structure of variables
%
%   AUTHOR: Jacob O. Wenegrat (jwenegrat@stanford.edu) 10/18/16
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thetas = 0:.005:2*pi; % Define angles to evaluate
r = (0.01:.01:1).*3.5*cr; % Define radial coordinates to evaluate

ubarmax = epsilon.*f.*cr; % Infer max velocity from inputs
jetwidth = cr;  % Assuming that the width of jet is set by eddy radius

% Gaussian SSH - Assuming a Gaussian SSH gives this vel profile:
velstruct = r./(jetwidth^2).*exp( - 1/2*(r/(jetwidth)).^2);
velstruct = velstruct/max(abs(velstruct));
velstruct = ubarmax.*velstruct;
dudr = (1-r.^2./jetwidth.^2).*velstruct./r; % Radial derivative of vel

disp(['Maximum \zeta_s = ', num2str(max(abs(dudr)))]);
disp(['Maximum \zeta_c = ', num2str(max(abs(velstruct./r)))]);
% Gaussian Vel - Uncomment to assume a Gaussian vel profile
% velstruct = ubarmax.*exp(- ((r-cr)/jetwidth).^2);
% dudr = -(r-cr).*2./(jetwidth).^2.*velstruct;

% PreAllocating variables
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
% maxom = 0;

% Iterate through radial position
for i=1:length(r)
   omega = velstruct(i)./r(i); % Curvature vorticity
%    if (abs(omega)>maxom); maxom = abs(omega); end
   zeta =  dudr(i) + omega; % Total vorticity
   
   % Following Equations 12-15 in Wenegrat and Thomas 2016
   denom = (f+2*omega).*(f+zeta)-omega.^2;
   
    % Note that these definitions are equivalent to the version given below. 
%    Mx(i,:) = (zeta - 2*omega)./denom.*sin(thetas).*cos(thetas).*tau;
%    My(i,:) = -(f + omega + 2.*omega.*sin(thetas).^2 + zeta.*cos(thetas).^2)./denom.*tau;
   
   Mr(i,:) = - (f + 3*omega)./denom.*sin(thetas).*tau;
   Mth(i,:) = -(f + omega+zeta)./denom.*cos(thetas).*tau;
 
   %Transport Equations
   Mx(i,:) = Mr(i,:).*cos(thetas) - Mth(i,:).*sin(thetas);
   My(i,:) = Mr(i,:).*sin(thetas) + Mth(i,:).*cos(thetas);
   
   Myc(i,:) = -tau./(f+zeta).*ones(size(thetas)); % Classic NL Ekman transport
   
   % Also doing the order epsilon expansions
   Myeps(i,:) = -(f - (zeta-omega).*sin(thetas).^2 - omega.*cos(thetas).^2).*tau./f.^2;
   Mxeps(i,:) = (zeta-2*omega).*sin(thetas).*cos(thetas).*tau; % Over f?
   
%    Mvclass(i,:) = - tau./(f+zeta); % Classic NL transport relation
   Mvlin(i,:) = -tau./f; % Classic Ekman transport
   
   [x(i,:), y(i,:)] = pol2cart(thetas, r(i)); % Convert to cartesian coordinates
   uc(i,:) = velstruct(i).*(-sin(thetas)); % Find the balanced flow in cart
   vc(i,:) = velstruct(i).*cos(thetas);
   umag(i,:) = velstruct(i).*ones(size(thetas)); %Not quite right, but works fine for these purposes
   omegaf(i,:) = omega; % Save the Omega field
   zetaf(i,:) = zeta;   % Save the zeta field
end

nr = length(r); nth = length(thetas);

xvec = reshape(x, nr*nth, 1); 
yvec = reshape(y, nr*nth, 1);
Mxvec = reshape(Mx, nr*nth, 1);
Myvec = reshape(My, nr*nth, 1);

x = (-1:.02:1).*r(end); y = (-1:.02:1).*r(end); % New coordinates for regridding
deltax = x(2)-x(1); deltay=y(2)-y(1);
[X, Y] = meshgrid(x, y);

% Regrid Ekman Transports
MUBAR = griddata(xvec, yvec, Mxvec, X, Y);
MVBAR = griddata(xvec, yvec, Myvec, X, Y);
[MUx, MUy] = gradient(MUBAR, deltax); % Tested this order, correct, 10/18/16
[MVx, MVy] = gradient(MVBAR, deltay);
W = MUx+MVy; % Ekman vertical vel

Mycvec = reshape(Myc, nr*nth, 1); %Classic NL Ek transport
MVBARC = griddata(xvec, yvec, Mycvec, X, Y);
[MVx, MVy] = gradient(MVBARC, deltay);
Wclassic = MVy; % Classic NL Pumping Vel

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
% mvclassvec = reshape(Mvclass,nr*nth,1);
% MVCLASS = griddata(xvec, yvec, mvclassvec, X, Y);
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

% Put output variables in struct
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
% out.MVCLASS= MVCLASS; %Not sure about this calculation, but not using it.
out.MVLIN = MVLIN;
end
