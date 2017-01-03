%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIRCULAR EKMAN TRANSPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Generates Figure 3 for NL Ekman Manuscript
% Following Analytical solution of Section 3
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% PARAMETERS
epsilon = 0.25;
% epsilon = 0.0266;
% epsilon = 0.029;
tau = .1/1035;
f=1e-4;
% f = 2*2*pi./86400*sind(30);
cr = 12e3; % Radius of Eddy
out = calcCircWGradientWind(-epsilon, tau, f,cr); % Calculate (epsilon < 0 is anticyclone)
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic; UMAG= out.UMAG;

%% PLOT ANTICYCLONE SIDE
figure
nc = 20;24;16; % Number of contours in plot
xl = [-2.5 2.5];
gap = [.03 .01]; margh = .1; margw = .1;
eps = max(max(abs(out.ZETA./f)));
eps = epsilon;
for i=1:4
    switch i
        case 1
            var = MUBAR./(eps*tau./f);
            conts =  linspace(-1, 1, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_x}{\epsilon\tau_o/(\rho f)}$';
            pos  = [3 5];
        case 2
            var = MVBAR./(tau./f);
            conts = linspace(-1 - eps, -1 + eps, nc);
            clim = conts([1 end]);
            conts = [-2 conts 0];
            lstr = '$\frac{M_y}{\tau_o/(\rho f)}$';
            pos = [7 9];
        case 3 
            var = W./(eps.*tau./(f*cr));
            conts = linspace(-0.25, 0.25, nc); %XX_CONFIRM THIS AGAIN
            conts = linspace(-10 , 10, nc);
            clim = conts([1 end]);
            lstr = '$\frac{w_e}{\epsilon\tau_o/(\rho fR)}$';
            pos = [11 13];
        case 4
            var = (W - Wclassic)./(eps.^2*tau./(f*cr));
            clim = conts([1 end]);
%             conts = [-0.3 linspace(-0.25, 0.25, nc) 0.3];
            lstr = '$\frac{w''}{\epsilon^2\tau_o/(\rho fR)}$';
            pos = [15 17];
%             conts = wcontsclass;
    end
subtightplot(9,2,pos, gap, margh, margw);
    [c, h] = contourf(X./cr, Y./cr, var, conts); shading interp
    % set(h, 'edgecolor', 'none');
    axis square
    % colorbar;
    set(gca, 'xlim', xl, 'ylim', xl);
    % grid on
    ylabel(lstr, 'Rotation', 0, 'Interpreter', 'Latex');
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
end
set(gcf,'Color', 'w', 'Position', [384   130   613   843])


%% VELOCITY STRUCTURE PLOT
ind = find(out.y == 0);

uvels = out.UMAG(:, ind);
uvels(1:ind) = -uvels(1:ind);
uvels(ind) = 0;
om = out.OMEGA(:,ind);
zetas = out.ZETA(:,ind) - om; % Shear vorticity

subtightplot(9, 2, 1, gap, margh, margw);
    plot(X./cr, uvels./(epsilon.*f.*cr), 'LineWidth', 2);
    hold on
        plot(X./cr, om./(f*epsilon*2), '--');
        plot(X./cr, zetas./(f*epsilon*2), '--');
    hold off
    % axis equal
    set(gca, 'XTick', -2:1:2, 'XLim', xl);
        % set(gca, 'ylim', [0 .4])
        % axis square
    grid on
    cb = colorbar;
    set(cb, 'visible', 'off');
    pos = get(gca, 'position');
    % set(gca, 'Position', [1.18 pos(2) 1.375 pos(4)]);
    set(gca, 'Position', [.191 pos(2) 0.212 pos(4)]);
    title('Anticyclone');
    set(gca, 'FontSize', 12);

%% CYLCONE
out = calcCircWGradientWind(epsilon, tau, f, cr);  % Recalculate velocity field
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic;
%%
for i=1:4
    switch i
        case 1
            var = MUBAR./(eps*tau./f);
            conts =  linspace(-1, 1, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_x}{\epsilon\tau_o/(\rho f)}$';
            pos  = [4 6];
        case 2
            var = MVBAR./(tau./f);
            conts = linspace(-1-eps, -1+eps, nc);
            clim = conts([1 end]);
            conts = [-2 conts 0];
            lstr = '$\frac{M_y}{\tau_o/(\rho f)}$';
            pos = [8 10];
        case 3 
            var = W./(eps.*tau./(f*cr));
            conts = linspace(-0.25, 0.25, nc); %XX_CONFIRM THIS AGAIN
            conts = linspace(-3 , 3, nc);
            clim = conts([1 end]);
            lstr = '$\frac{w_e}{\epsilon\tau_o/(\rho fR)}$';
            pos = [12 14];
        case 4
            var = (W - Wclassic)./(eps.^2*tau./(f*cr));
            clim = conts([1 end]);
%             conts = [-0.3 linspace(-0.25, 0.25, nc) 0.3];
            lstr = '$\frac{w''}{\epsilon^2\tau_o/(\rho fR)}$';
            pos = [16 18];
%             conts = wcontsclass;
    end
subtightplot(9,2,pos, gap, margh, margw);
    [c, h] = contourf(X./cr, Y./cr, var, conts); shading interp
    % set(h, 'edgecolor', 'none');
    axis square
    cb = colorbar;
    if i==2
       set(cb, 'Ticks', [clim(1) -1 clim(end)]);
    else
        set(cb, 'Ticks', [clim(1) 0 clim(end)]);
    end
    set(gca, 'xlim', xl, 'ylim', xl);
    % grid on
    % ylabel(lstr, 'Rotation', 0, 'Interpreter', 'Latex');
    % hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
    hold on; 
        plot(cos(0:.1:2*pi), sin(0:.1:2*pi), '--k'); 
        plot(2*cos(0:.1:2*pi), 2.*sin(0:.1:2*pi), '--k'); 
    hold off
    grid on
    set(gca, 'clim', clim);
    set(gca, 'XTick', -2:1:2, 'YTick', -2:1:2);
    set(gca, 'FontSize', 12);
end
colormap(flipud(othercolor('RdBu11')));
set(gcf,'Color', 'w', 'Position', [384   130   613   843])

%% CYCLONE VELOCITY PLOT

subtightplot(9,2, 2, gap, margh, margw);
    ind = find(out.y == 0);
    uvels = out.UMAG(:, ind);
    uvels(1:ind) = -uvels(1:ind);
    uvels(ind) = 0;

    om = out.OMEGA(:,ind);
    zetas = out.ZETA(:,ind) - om;
    plot(X./cr, uvels./(epsilon.*f.*cr), 'LineWidth', 2);
    hold on
        plot(X./cr, om./(f*epsilon*2), '--');
        plot(X./cr, zetas./(f*epsilon*2), '--');
    hold off
    % axis equal
    set(gca, 'XTick', -2:1:2, 'XLim', xl);
    % set(gca, 'ylim', [0 .4])
    % axis square
    grid on
    cb = colorbar;
    set(cb, 'visible', 'off');
    pos = get(gca, 'position');
    set(gca, 'Position', [.5408 pos(2) 0.2150 pos(4)]);
    title('Cyclone');
    set(gca, 'FontSize', 12);
%     legend('$\bar{u}/(\epsilon f R)$', '$\Omega/(2\epsilon f)$', '$(\zeta-\Omega)/(2\epsilon f)$', 'Location', 'EastOutside');
    leg = legend('$\frac{\bar{u}}{\epsilon f R}$', '$\frac{\Omega}{2\epsilon f}$', '$\frac{\zeta-\Omega}{2\epsilon f}$', 'Location', 'EastOutside');
    set(leg, 'FontSize', 12);
    %% EXPORTING - TO EXPORT USE THIS COMMAND:
% First fix Legend location, and ylabel of bottom row.
% export_fig('CircularWNew.eps', '-eps', '-opengl', '-q100', '-p01')