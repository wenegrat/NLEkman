%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GaussBumpWithVorticity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Makes Figure 6 for NL Ekman manuscript
% Calculates vorticity terms for a front undergoing a Gaussian meander
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Define parameters
deltax = 250;
x=0:deltax:2500e3;
f = 1e-4;
% ubar = .25;
% xc = floor(length(x)/2); % Find center point
% xc = x(xc);
xc = 499750+5e5; 
amp = 12e3; %Amplitude of bump
% amp = 6e3;

width = amp/1.5; % Width of bump
y = amp*exp( - ((x-xc)/width).^2); % Define path


% % Alternate is Cauchy distribution
% gamma = width;
% y =  (1./(pi*gamma* ( 1+ ( (x-xc)./gamma).^2)));
% y= y./max(y);
% y = amp.*y;

plot(x, y); % Basic frontal structure plot

% Ramp function for wind stress
xc = 420e3;
rampwidth = 120e3;
facTau = 1/2.*(1+tanh( (x-xc)./rampwidth));
% facTau = 1;
tau = .1.*facTau.*ones(size(x))/1035;
tau = tau-tau(1);

epsilon = .5; % Define epsilon
ubarc = .5; % Max balanced vel
L = ubarc./(f*epsilon); % Infer L scale
deltay = 500; 

offsets = -amp/2:deltay:2*amp; % Positions to calculate solution at
%Gaussian Vel profile
ubart = ubarc.*exp(-(offsets./L).^2/2);

%  % Linear ubart profile
% hw = 13;
% % ubart = NaN(size(offsets));
% dudnt = ubart;
% ubart(1:hw) = ubarc.*(1+(offsets(1:hw)./offsets(1)).*(ubart(1)-ubarc)./ubarc);
% ubart(hw:end) = ubarc.*(1+((offsets(hw:end)-offsets(hw))./offsets(end)).*(ubart(end)-ubarc)./ubarc);
% dudnt(1:hw) = (1./offsets(1)).*(ubart(1)-ubarc);
% dudnt(hw:end) = (1./offsets(end)).*(ubart(end)-ubarc);
% dudnt(hw) = 0;

% Pre -Allocate
u = NaN(length(offsets), length(x));
v = u; l = u; zeta = u; k = u; taus = u; taun = u; positions=u;ux=u; vx=u; 
svec = u; ubarall = u; omega = zeta;

% Solve along each offset curve
for i=1:length(offsets)
    disp([num2str(i),'/', num2str(length(offsets))]);
    off=offsets(i);
    pos = makeOffsetCurve(y, off, x);
    if isfinite(pos); 
        ubar = ubart(i);
        dudn = -ubar.*offsets(i)./(L.^2); % GAUSSIAN PROFILE
    %     dudn = -dudnt(i); % LINEAR PROFILE
     
        % Solution of ODE takes place here:
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        out = SolveFrontEkman(x, pos, ubar, -dudn, tau, f);
        
        % Save output variables
        u(i,:) = out.u;
        ux(i,:) = out.ux;
        v(i,:) = out.v;
        vx(i,:) = out.vy;
        l(i,:) = out.l;
        zeta(i,:) = out.zeta;
        omega(i,:) = out.omega;
        k(i,:) = out.omega./ubar;
        taus(i,:) = out.taus;
        taun(i,:) = out.taun;
        positions(i,:) = pos;
        svec(i,:) = out.tangent;
        ubarall(i,:) = ubar.*ones(size(x));
    end
end

%% INTERPOLATE TO GRID
xi = repmat(x, [1, length(offsets)]).';
noff = length(offsets); nx = length(x);
% lvec = reshape(l, noff*nx,1);
yi = reshape(positions.', noff*nx,1);
uvec = reshape(ux.', noff*nx, 1);
vvec = reshape(vx.', noff*nx, 1);
uevec = reshape(u.', noff.*nx, 1);
vevec = reshape(v.', noff*nx,1);
svecvec = reshape(svec.', noff*nx, 1);
kvec = reshape(k.', noff*nx,1);
ubarvec = reshape(ubarall.', noff*nx, 1);
zetavec = reshape(zeta.', noff*nx, 1);
omegavec = reshape(omega.', noff*nx, 1);

[X, Y] = meshgrid(linspace(x(1), x(end), 4000), offsets);
% [X, Y] = meshgrid(x, offsets);

deltax = X(1,2) - X(1,1);
deltay = Y(2,1) - Y(1,1); %

masky = isfinite(yi);
U = griddata(xi(masky), yi(masky), uvec(masky), X, Y);
V = griddata(xi(masky), yi(masky), vvec(masky), X, Y);
SVECX = griddata(xi(masky),yi(masky), svecvec(masky), X, Y);
K = griddata(xi(masky), yi(masky), kvec(masky), X, Y);
UBAR = griddata(xi(masky), yi(masky), ubarvec(masky), X, Y);
ZETABAR = griddata(xi(masky), yi(masky), zetavec(masky), X, Y);
OMEGABAR = griddata(xi(masky), yi(masky), omegavec(masky), X, Y);
Ue = griddata(xi(masky), yi(masky), uevec(masky), X, Y);
Ve = griddata(xi(masky), yi(masky), vevec(masky), X, Y);

mask = ones(size(Y)); 
ytemp = interp1(x, y, X(1,:));

% Create a mask for positions with offsets too large for the curvature
for i=1:length(ytemp);
    mask(Y(:,i)> ytemp(i)+amp, i) = NaN;
    mask(Y(:,i)< ytemp(i)+offsets(1), i) = NaN;
end
[Ux, Uy] = gradient(U, deltax, deltay);
[Vx, Vy] = gradient(V, deltax, deltay);

Ux = Ux.*mask; Uy = Uy.*mask; Vx=Vx.*mask; Vy=Vy.*mask;
WE = Ux + Vy; ZETABAR  = ZETABAR.*mask;

VortCart = Vx-Uy;


dvds = real(SVECX).*Vx + imag(SVECX).*Vy;
dudn = -imag(SVECX).*Ux + real(SVECX).*Uy;
uk = (real(SVECX).*U + imag(SVECX).*V).*K;
VortNat = dvds - dudn + uk;

%% CALCULATE TERMS IN BUDGET
% Basic approach is to calculate terms in cartesian coordinates and later
% interpolate along a path.

%Along front advection of Ekman vorticity
[ZetaX, ZetaY] = gradient(VortCart, deltax, deltay);
ADV = UBAR.*( real(SVECX).*ZetaX + imag(SVECX).*ZetaY);

STRETCH = -(f+ZETABAR).*WE; % This is correct sign (as we are considering the vertical integral from 0 to -z).

[ZBarX, ZBarY] = gradient(ZETABAR, deltax, deltay);
GRAD = -Ue.*( real(SVECX).*ZBarX + imag(SVECX).*ZBarY) - Ve.*(-imag(SVECX).*ZBarX + real(SVECX).*ZBarY);
GRADS = -Ue.*( real(SVECX).*ZBarX + imag(SVECX).*ZBarY);
GRADN = - Ve.*(-imag(SVECX).*ZBarX + real(SVECX).*ZBarY);

%Checking just the curvature vort component.
[ZBarOX, ZBarOY] = gradient(OMEGABAR, deltax, deltay);
GRADNK = - Ve.*(-imag(SVECX).*ZBarOX + real(SVECX).*ZBarOY);

%% Make Figure 6
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot Parameters
xnorm = 2*pi*ubarc./(f);
lxpos = .125;
fs = 12;
% xl = [4.24e5 8e5]+5e5;
xl = [4.25e5 8e5]+5e5; % Plot x limits

index = find(x>xl(1), 1)-1;
indexe = find(x>xl(2), 1);
xl = [0 6]; % Convert to normalized units
% xl = [0 30];
xn = x./xnorm;
xn = xn - xn(index); % start at 0 to make plots nice
gap = [.02 .05]; margh = .1; margw=.1;

% Pick path to plot
ind = 13; %11 works well too, 13 is center of jet, so probably simplest to do this.

% Interpolate to along meander
ADVI = interp2(X, Y, ADV, x, positions(ind,:));
STRETCHI = interp2(X, Y, STRETCH, x, positions(ind,:));
GRADI = interp2(X,Y, GRAD, x, positions(ind,:));

% Make Plots
figure
subtightplot(3,1,1,gap, margh, margw)
    rp = 6; % Only plot some of the vectors
    plot(xn, positions(ind,:)./xnorm, 'LineWidth', 2);
    hold on
        set(gca, 'ColorOrderIndex', 1);
        q = quiver(xn(1:rp:end), positions(ind,1:rp:end)./xnorm, ux(ind,1:rp:end), vx(ind,1:rp:end));
        set(q, 'ShowArrowHead', 'off');
    hold off
    axis equal % necessary for stick vectors
    set(gca, 'xlim', xl);
    grid on
    % title(num2str(Y(ind,1)));
    % axis equal
    set(gca, 'ylim', [-2 .5]);
    set(gca, 'XTickLabel', []);
    ylabel('$\hat{y}$', 'Interpreter', 'Latex');
    t = text(lxpos, 0.25, 'Ekman Transport');
    set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
    
    set(gca, 'FontSize', fs);

% vortnorm = taumag./(f.*L).*ubarc;
vortnorm = ubarc.*f.*20/L;
subtightplot(3,1,2,gap, margh, margw)
    plot(xn, ADVI./vortnorm, 'LineWidth', 2);
    hold on
    plot(xn, STRETCHI./vortnorm, 'LineWidth', 2);
    plot(xn, GRADI./vortnorm, 'LineWidth', 2);
    hold off
    set(gca, 'xlim', xl);
    set(gca, 'YTick', -1.5:0.5:1.5);
    grid on
    % set(gca,'ylim', [-10 20]);
    set(gca, 'ylim', [-1 1.5]);
    legend('ADV', 'STRETCH', 'GRAD', 'Location', 'NorthEast');
    set(gca, 'XTickLabel', []);
    t = text(lxpos, 1.25, 'Vorticity Budget');
    set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
    set(gca, 'FontSize', fs);

vortnorm = vortnorm.*L;
subtightplot(3,1,3,gap, margh, margw)
plot(xn(index:indexe), cumtrapz(x(index:indexe), ADVI(index:indexe))./vortnorm, 'LineWidth', 2);
hold on
    plot(xn(index:indexe), cumtrapz(x(index:indexe), STRETCHI(index:indexe))./vortnorm, 'LineWidth', 2);
    plot(xn(index:indexe), cumtrapz(x(index:indexe), GRADI(index:indexe)./vortnorm), 'LineWidth', 2);

    % This can be used as confirmation that budget closes.
    % plot(xn(index:indexe), cumtrapz(x(index:indexe), ADVI(index:indexe) - ...
    %     STRETCHI(index:indexe)-GRADI(index:indexe))./vortnorm, 'LineWidth', 2, 'LineStyle', '--');
hold off
set(gca, 'xlim', xl);
set(gca, 'ylim', [-1.5 1.5]);
set(gca, 'YTick', -1.5:0.5:1.5);

grid on
xlabel('$\hat{x}$', 'Interpreter', 'Latex');
t = text(lxpos, 1.23, 'Integrated Budget');
set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
set(gca, 'FontSize', fs);
set(gcf,'Color', 'w', 'Position', [650   253   521   706])

% TO EXPORT, FIX LEGEND POSITION THEN
% export_fig('VorticityBudgetNew.eps', '-eps', '-painters', '-q101', '-p01')

%%
% %% Plot Ekman Vorticity Along the path choosen for the final line plot
% VI = interp2(X, Y, VortCart, x, positions(ind,:));
% figure 
% plot(x./xfact, VI);
%% CONFIRM VORTICITY
% % Show that the calculation of vorticity in nat coordinates is equivalent
% % with cartesian.
% cl = [-1 5].*1e-4;
% figure
% subplot(2,1,1)
% pcolor(X,Y, VortCart); shading interp; colorbar;
% set(gca, 'clim', cl);
% %     set(gca, 'xlim', [4.4e5 5.6e5]);
% hold on
%     plot(x, positions(ind,:))
% hold off
% subplot(2,1,2)
% pcolor(X,Y, VortNat); shading interp; colorbar;
% set(gca, 'clim', cl);
% %     set(gca, 'xlim', [4.4e5 5.6e5]);
% hold on
%     plot(x, positions(ind,:))
% hold off

%% CONTOURS
% % Contour Plots of the terms in the vort budget
% % Useful for testing closure of budget
% 
% cl = [-1 1].*5e-8;
% 
% mask = ones(size(Y)); 
% for i=1:length(x)
% %     mask(Y(:,i)> y(i)+offsets(end)/4, i) = NaN;
% %     mask(Y(:,i)< y(i)+offsets(1)/1, i) = NaN;
% end
% conts = linspace(cl(1), cl(end), 10);
% figure
% for i=1:5
%     switch i
%         case 1
%             VAR = VortCart.*1e-4;
%         case 2
%             VAR = ADV;
%         case 3
%             VAR = STRETCH;
%         case 4
%             VAR = GRAD;
%         case 5
%             VAR = ADV-STRETCH-GRAD;
%     end
%    subplot(5,1,i) 
%     [c, h] = contourf(X, Y, VAR.*mask, conts); shading interp
%     set(h, 'edgecolor','none')
%     colorbar;
%     set(gca, 'clim', cl);
% 
% %     set(gca, 'xlim', [4.8e5 7.8e5]);
% %     set(gca, 'ylim', [offsets(1) offsets(end)]);
%         hold on
%     ind = 11;
%     plot(x, positions(ind,:), 'k'); 
% %     quiver(x, positions(ind,:), ux(ind,:), vx(ind,:), .25);
%     hold off
% %     axis equal
% end



