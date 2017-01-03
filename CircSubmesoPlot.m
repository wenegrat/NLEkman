%% PARAMETERS
epsilon = 0.25;
% epsilon = 0.3;
% epsilon = 0.0266;
% epsilon = 0.029;
tau = .1/1035;
% f=1e-4;
% f = 2*2*pi./86400*sind(30);
f= 1e-4;
cr = 12e3; % Radius of Eddy
out = calcCircWGradientWind(-epsilon, tau, f,cr); % Calculate (epsilon < 0 is anticyclone)
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic; UMAG= out.UMAG;


%% PLOT ANTICYCLONE SIDE
denom =  (f+out.ZETA).*(f+2*out.OMEGA) - out.OMEGA.^2;

figure
pcolor(X./cr,Y./cr, W./(epsilon*tau./(f*cr))); shading interp
pcolor(X./cr,Y./cr, denom); shading interp
cl = get(gca, 'clim');
hold on
contour(X./cr, Y./cr, denom, [0 0],'k');
% contour(X./cr,Y./cr,out.ZETA./f, [-1 1]./2, '--k', 'LineWidth', 2);
% contour(X./cr,Y./cr,abs(out.ZETA-out.OMEGA)./f, [0.5 0.5], 'r');
% contour(X./cr,Y./cr,out.ZETA./f, [-1 1], 'k');
% contour(X./cr,Y./cr,(out.ZETA + out.OMEGA)./f, [-1 -1], 'b', 'LineWidth', 2);
% contour(X./cr,Y./cr,out.OMEGA./f, [-1 1]/3, 'LineWidth', 2)
hold off
set(gca, 'clim',cl);
%%
eps = max(max(abs(out.ZETA./f)))
eps = 0.6
umax = max(max(abs(out.UMAG)))
% eps = epsilon
figure
nc = 20;24;16; % Number of contours in plot
xl = [-2.5 2.5];
gap = [.125 .025]; margh = .1; margw = .1;
fs=20;
for i=1:4
    switch i
        case 1
            var = MUBAR./(eps*tau./f);
            conts =  linspace(-1, 1, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_x}{\epsilon\tau_o/(\rho f)}$';
            ts = 'A)';
            pos  = 1;
        case 2
            var = MVBAR./(tau./f);
            conts = linspace(-1-eps, -1 + eps, nc);
            clim = conts([1 end]);
            conts = [-100 conts 100];
            lstr = '$\frac{M_y}{\tau_o/(\rho f)}$';
            pos = 3;
                        ts = 'B)';

        case 3 
            var = W./(eps.*tau./(f*cr));
            conts = linspace(-0.25, 0.25, nc); %XX_CONFIRM THIS AGAIN
            conts = linspace(-3, 3, nc);
            clim = conts([1 end]);
%             conts = [-100  conts 100]
            lstr = '$\frac{w_e}{\epsilon\tau_o/(\rho fR)}$';
            pos = 2;
                        ts = 'C)';

        case 4
            var = (W - Wclassic)./(eps.^2*tau./(f*cr));
            clim = conts([1 end]);
%             conts = [-0.3 linspace(-0.25, 0.25, nc) 0.3];
            lstr = '$\frac{w''}{\epsilon^2\tau_o/(\rho fR)}$';
            pos = 4;
                        ts = 'D)';

%             conts = wcontsclass;
    end
subtightplot(2,2,pos, gap, margh, margw);
    [c, h] = contourf(X./cr, Y./cr, var, conts); shading interp
    % set(h, 'edgecolor', 'none');
    axis square
    cb = colorbar;
    set(cb, 'TickLabelInterpreter', 'Latex')
    if i==2
       set(cb, 'Ticks', [clim(1) -1 clim(end)]);
    else
        set(cb, 'Ticks', [clim(1) 0 clim(end)]);
    end
    set(gca, 'xlim', xl, 'ylim', xl);
    % grid on
    set(get(gca, 'ylabel'), 'Position', [-4.5260   -0.5504         0]);
    % hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
    hold on; 
        plot(cos(0:.1:2*pi), sin(0:.1:2*pi), '--k'); 
        plot(2.*cos(0:.1:2*pi), 2.*sin(0:.1:2*pi), '--k'); 
    hold off
    grid on
    set(gca, 'clim', clim);
    set(gca, 'XTick', -2:1:2, 'YTick', -2:1:2);
    set(gca, 'FontSize', 12);
    set(get(gca, 'ylabel'), 'FontSize', 18);
        title(lstr, 'Rotation', 0, 'Interpreter', 'Latex', 'FontSize', fs);
    t = text(-2.2, 2, ts);
        set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
end
% set(gcf,'Color', 'w', 'Position', [384   130   613   843])
colormap(flipud(othercolor('RdBu11')));
set(gcf,'Color', 'w', 'Position', [   392     3   746   545])

% export_fig('SubmesoscaleVortex.eps', '-eps', '-opengl', '-q101', '-p01')
