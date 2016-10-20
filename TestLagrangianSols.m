%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% TestLagrangianSols.m
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Try solving the Lagrangian EOM, and comparing to our streamwise solution.
%
% Provides a basic confirmation of the oscillatory behavior we see.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% PARAMETERS
f=1e-4;
epsilon = 0.25;
ubarmax = 1;
L = ubarmax./(f*epsilon);

deltax= 500;
x = 0:deltax:(120e3*30);
amp = 4e3;
rampwidth = 120e3.*4;
xc = 120e3*8;
facA = 1/2.*(1+tanh( (x-xc)./rampwidth));
ypos = facA.*amp.*sin(2*pi.*x./120e3);

deltay = 500;
ys = -3*amp:deltay:3*amp; %Offset positions
utot = ubarmax.*exp(-(ys./L).^2/2);
dudntot = -ys./L.^2.*utot;

% Pre-Allocate
ytotsl = NaN(length(ys), length(x));
ux = ytotsl; uy = ytotsl; ktot = ytotsl;
% Iterate over offsets
for i=1:length(ys)
   y = makeOffsetCurve(ypos, ys(i), x);
    du = gradient(x, deltax); %Temporary velocity gradients
    dv = gradient(y, deltax); %
    dv(1:end-1) = (y(2:end)-y(1:end-1))./(deltax);
    dv(end)=dv(1);
    vels = du+1i*dv;  %Tangent Vectors at each spot.
    frntvec = vels./abs(vels);
    
        if ys(i) == 0;
            velt = vels;
        end
        
    ytotsl(i,:) = y;
    
    % Project components of imposed balanced flow on cart coordinates
    ux(i,:) = utot(i).*real(frntvec);
    uy(i,:) = utot(i).*imag(frntvec);
    
    % Calculate curvature
    dx  = gradient(x, deltax);
    ddx = gradient(dx, deltax);
    dy = gradient(y, deltax);
    ddy = gradient(dy, deltax);
    dy(1:end-1) = (y(2:end)-y(1:end-1))./(deltax);
    dy(end) = dy(1);
    ddy(1:end-1) = (dy(2:end)-dy(1:end-1))./(deltax);
    ddy(end) =ddy(1);
    num   = dx .* ddy - ddx .* dy;
    denom = dx .* dx + dy .* dy;
    denom = sqrt(denom);
    denom = denom.* denom.* denom;
    k = num ./ denom;
    k(denom < 0) = NaN;
    ktot(i,:) = k;
    
end

[ny, nx] = size(ytotsl);
yvec = reshape(ytotsl.', ny*nx, 1);
xvec = repmat(x, [1 ny]).';
uxvec = reshape(ux.', ny*nx,1);
uyvec = reshape(uy.', ny*nx, 1);
kvec = reshape(ktot.', ny*nx, 1);

[X, Y] = meshgrid(x, ys);
U = griddata(xvec, yvec, uxvec, X, Y);
V = griddata(xvec, yvec, uyvec, X, Y);
K = griddata(xvec, yvec, kvec, X, Y);

[Ux, Uy] = gradient(U, deltax, deltay);
[Vx, Vy] = gradient(V, deltax, deltay);
%% DEFINE TAU FIELD

% tau is time-dependent, need to ramp up slowly to minimize oscillations
deltat = 0.25*3600;
t = 0:deltat:(24*3600*50);
[ny nx] = size(X);

%Note dividing by h for slab layer model, need to choose the parameters to satisfy epsilon_e<<1
tau = 1*.025./1035.*ones(length(t),ny, nx)./50; 

rampwidth = 24*3600*1;
tc = 24*3600*4;
facTau = 1/2.*(1+tanh( (t-tc)./rampwidth));
%         facTau = 1 + 0*facTau;
% facTau = ones(size(t));
tau = repmat(facTau.', [1 ny nx]).*tau; % Define tau uniform everywhere
tau = tau-tau(1,1,1);
%  tau = 0*tau;
%% DEFINE PRESSURE FIELD

% PY = f.*U + U.^2.*K;
% PX = -f.*V + V.^2.*K;

% Cyclogeostrophic Balance, note (PX, PY defined with sign convention such that they are the
% PGF)
PY = f.*U + U.*Vx + V.*Vy;
PX = -f.*V + U.*Ux + V.*Uy;

%% SOLVE LAGRANGIAN EOMS

initpos = Y(:,1);

xtots = NaN(length(initpos), length(t));
ytots = xtots; utots=xtots; vtots=xtots; ttots = xtots;
for j=1:length(initpos)
    guess = [0 U(j,1) initpos(j) V(j,1)];
    out = LagrangianIVP(t, guess, f, PX, PY, tau, X, Y); % xx- check this
    xtots(j,:) = out.x;
    ytots(j,:) = out.y;
    utots(j,:) = out.u;
    vtots(j,:) = out.v;
    ttots(j,:) = out.t;
end

% GRIDDATA
nt = length(t);
xvecs = reshape(xtots.', nt.*ny,1);
yvecs = reshape(ytots.', nt.*ny,1);
uvecs = reshape(utots.', nt.*ny,1);
vvecs = reshape(vtots.', nt.*ny,1);
tvecs = reshape(ttots.', nt.*ny, 1);

mask = isfinite(xvecs+yvecs+uvecs+vvecs);

UT = griddata(xvecs(mask), yvecs(mask), uvecs(mask), X, Y);
VT = griddata(xvecs(mask), yvecs(mask), vvecs(mask), X, Y);
TT = griddata(xvecs(mask), yvecs(mask), tvecs(mask), X, Y);

%% ALONG A PATH, SOLVE STREAMWISE PROBLEM

ind = 20; % offset position
uek = interp2(X, Y, UT-U, x, ytotsl(ind,:)); % Find the Lagrangian vals
vek = interp2(X, Y, VT-V, x, ytotsl(ind,:));
% utot = interp2(X, Y, abs(U+1i*V), x,ypos);

ttot = interp2(X,Y, TT, x, ytotsl(ind,:)); 

l = cumtrapz(x,abs(velt));

% Note there is some mixing of time/space dependence here
ti = interp1(t, squeeze(tau(1:end,1,1)), ttot); %Not qutie right, should use l?
% xpos = interp2(xtot
mask = isfinite(xtots(27,:));
% ti = interp1(xtots(27,mask), squeeze(tau(mask,1,1)), x); %Not qutie right, should use l?
% ti = interp1(ubarmax.*t, squeeze(tau(:, 1, 1)), x);
zetabase = -dudntot(ind);
out =    SolveFrontEkman(x, ytotsl(ind,:), ubarmax, zetabase,ti,  f);

xi = interp1(x, l, x);
% vyi = interp1(l, out.vy, x);
vyi = interp2(X, Y, griddata(x, ytotsl(ind,:), out.vy, X, Y),x, ytotsl(ind,:));

% Make Plot comparing solutions
figure
subplot(2,1,1);
plot(x./120e3, uek, 'LineWidth', 2);
hold on
plot(x./120e3, out.ux,'-', 'LineWidth', 1.5);
hold off
set(gca, 'ylim', [-1 1].*0.005);
title(['$\epsilon$ = ', num2str(epsilon), '  $\bar{u}$ = ', num2str(ubarmax), '  A = ', num2str(amp), '  $\zeta_s$ = ', num2str(zetabase)]);
legend('Lagrangian', 'Nonlinear Ekman Solution');
ylabel('u');
subplot(2,1,2);
plot(x./120e3, vek, 'LineWidth', 2);
hold on
plot(x./120e3, out.vy, 'LineWidth', 1.5);
hold off
set(gca, 'ylim', [-0.01 0]);
ylabel('v')
%% CONTOUR PLOTS SHOWING DOMAIN OF LAGRANGIAN SOLS
subplot(2,1,1);
pcolor(X./120e3, Y, UT-U); shading interp
hold on
plot(x./120e3, ypos); 
hold off
subplot(2,1,2);
pcolor(X./120e3, Y, VT-V); shading interp
hold on
plot(x./120e3, ypos); 
hold off
%% ZONAL MOM BALANCE
% Useful for understanding what is forcing the oscillations

dudt = gradient(interp2(X, Y, U, x, ytotsl(ind,:)), deltat);
cori = -f.*interp2(X, Y, V, x, ytotsl(ind,:));
pxi = interp2(X, Y, PX, x, ytotsl(ind,:));
tx = ones(size(pxi)).*tau(end,1,1);
figure
subplot(2,1,1)
plot(dudt);
hold on;
plot(cori);
plot(pxi, '--');
plot(tx);
% plot(dudt+cori, '--');
hold off
subplot(2,1,2)
plot(dudt);
hold on;
plot(pxi-cori);
hold off;