% A = Front_Solver(l, xfact, omega, zeta, f, ubar);
% ef = eig(A);
% 
% %Exponential growth rate
% nu = log(max(abs(ef(1)),abs(ef(2))));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
 % Should non-dimensionalize the problem so that I can vary just f/ubar and
 % epsk.
 
xfact = 120e3;
deltax = 1e3;
x = 1:deltax:120e3;
f = 1e-4;

zetabase = 0*.1*f.*ones(size(x)); zeros(size(x));
avecs = (0:0.1:24)*1e3;
kratio = .1:.01:6;
ubars  = f*xfact./(2*pi*kratio);
tic
nu = NaN(length(avecs), length(ubars));
ls = NaN(size(avecs));
kall = ls;
for i=1:length(avecs)
    disp(num2str(i))
    amp = avecs(i);
    y = amp*sin(2*pi.*x./xfact);

    %Determine curvature - Probably could do this analytically.
    dx  = gradient(x, deltax);
    ddx = gradient(dx, deltax);
    dy  = gradient(y, deltax);
    dy(1:end-1) = (y(2:end)-y(1:end-1))./(deltax);
    dy(end) = dy(1);
    ddy = gradient(dy, deltax);
    ddy(1:end-1) = (dy(2:end)-dy(1:end-1))./(deltax);
    ddy(end) =ddy(1);
    num   = dx .* ddy - ddx .* dy;
    denom = dx .* dx + dy .* dy;
    denom = sqrt(denom);
    denom = denom.* denom.* denom;
    k = num ./ denom;
    k(denom < 0) = NaN;
    kall(i) = max(k);
    
    vels = dx+1i*dy;  %Tangent Vectors at each spot.
    l = cumtrapz(x, abs(vels));
    ls(i) = l(end);
    parfor j=1:length(ubars)
        ubar = ubars(j);
    
        if (f/ubar <= 3*max(k))
            nu(i,j) = NaN;
        else
        omega = k*ubar;
        zeta = zetabase + omega;
%         A = Front_Solver(l, l(end), omega, zeta, f, ubar);
        A = Front_Solver(x, xfact,l, omega, zeta, f, ubar);
%         A = Front_Solver(linspace(0,2*pi, length(x)), 2*pi, omega, zeta, f, ubar);

        ef = eig(A);
        nu(i,j) = log(max(abs(ef(1)), abs(ef(2))));
        end
    end
end
toc

%%
k = f./ubars;
clear kn;
for j=1:length(ubars)
   kn(:,j) = k(j)./(2*pi./ls);
end
amat = repmat((avecs./ls).', [1 length(ubars)]);
amat = repmat((kall./(2*pi./ls)).', [1 length(ubars)]);

% amat(~isfinite(amat))=0;
knvec = reshape(kn, length(avecs).*length(ubars), 1);
amvec = reshape(amat, length(avecs).*length(ubars),1);
nuvec = reshape(nu, length(avecs).*length(ubars),1);
[X,Y ] = meshgrid(k./(2*pi./xfact), kall./(2*pi./xfact));
NU = griddata(knvec, amvec, nuvec, X, Y);
% %%
% figure
% pcolor(X, Y, NU); shading interp
% hold on
% cl = get(gca, 'Clim');
% contour(X, Y, NU, [.1 .1], 'k')
% set(gca, 'clim', cl);
% hold off
% xlabel('$(f/\bar{u})/(2\pi/\lambda_l)$', 'Interpreter', 'Latex')
% ylabel('$A/\lambda_l$', 'Interpreter', 'Latex');
% colormap(othercolor('Reds9'));
% colorbar
% set(gcf, 'Color', 'w');
% grid on
%%
x=0:1e3:1000e3;
ubar = .25;
xc = floor(length(x)/2);
xc = x(xc);
amp = 15e3;
width = amp/2;
y = amp*exp( - ((x-xc)/width).^2);

plot(x, y);
xc = 120e3;
rampwidth=80e3;
facTau = 1/2.*(1+tanh( (x-xc)./rampwidth));
% facTau = 1;
% tau = facTau.*taumag.*ones(size(x));

tau = .1.*facTau.*ones(size(x));
tau = tau-tau(1);
zetabase = .1*f;
out = SolveFrontEkman(x, y, ubar, zetabase, tau, f);
subplot(4,1,1)
plot(x, tau);
subplot(4,1,2:4)
quiver(x, y, out.ux, out.vy);
hold on; plot(x, y, 'LineWidth', 2); hold off

nrm = 2*pi./(f./ubar);

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
cb = colorbar;
set(cb, 'visible', 'off')
% axis equal
hold on
% quiver(.5, .75, tau(end), 0)
hold off
grid on
xlabel('$x/\lambda_n$', 'Interpreter', 'Latex');
ylabel('$y/\lambda_n$', 'Interpreter', 'Latex');
set(gca, 'FontSize', 20);
%%
% feff = mean(sqrt((f+2*omega).*(f+zeta)));
feff = f;
% figure
subplot(2,1,2)
pcolor(f./ubars./(2*pi./xfact), avecs./xfact, nu); shading interp;
pcolor((feff)./ubars./(2*pi./xfact), kall./(2*pi./xfact), nu); shading interp;

% pcolor(f./ubars, avecs, nu); shading interp;

hold on
cl = get(gca, 'Clim');

% contour((f./ubars./(2*pi./xfact)), avecs./xfact, nu, [.1 .1], 'k')
% contour((f./ubars./(2*pi./xfact)), kall./(2*pi./xfact), nu, [.1 .1]/2, 'k')

set(gca, 'clim', cl);
hold off
xlabel('$(f/\bar{u})/(2\pi/\lambda_l)$', 'Interpreter', 'Latex')
% ylabel('$A/\lambda_x$', 'Interpreter', 'Latex');
ylabel('$(\Omega_{max}./\bar{u})/(2\pi/\lambda_l)$', 'Interpreter', 'Latex');

colormap(othercolor('Reds9'));
cb = colorbar
set(get(cb, 'ylabel'), 'String', 'Growth Rate (\lambda_l^{-1})');
set(gcf, 'Color', 'w');
set(gca, 'FontSize', 20);
grid on
set(gca, 'clim', [0 0.25])