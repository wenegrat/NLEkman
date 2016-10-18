%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CIRCULAR EKMAN TRANSPORT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
epsilon = 0.05;
tau = .1/1035;
f=1e-4;
cr = 75e3;
out = calcCircW(-epsilon, tau, f,cr);
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic; UMAG= out.UMAG;
%%
% vlim = [-1.25 -.85];
% wlim = [-1 1]*51e-5;
% wlim = [-1 1].*1e-5;


%%
figure
nc = 30;20;24;16;

xl = [-2.5 2.5];

gap = [.03 .01]; margh = .1; margw = .1;

for i=1:4
    switch i
        case 1
            var = MUBAR./(epsilon*tau./f);
            conts =  linspace(-1.2, 1.2, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_x}{\epsilon\tau_o/(\rho f)}$';
            pos  = [3 5];
        case 2
            var = MVBAR./(tau./f);
            conts = linspace(-1.2, -.8, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_y}{\tau_o/(\rho f)}$';
            pos = [7 9];
        case 3 
            var = W./(tau./(f*cr));
            conts = linspace(-0.25, 0.25, nc);
            clim = conts([1 end]);
            lstr = '$\frac{w_e}{\tau_o/(\rho fR)}$';
            pos = [11 13]
        case 4
            var = (W - Wclassic)./(2*epsilon*tau./(f*cr));
            clim = conts([1 end]);
            conts = [-0.3 linspace(-0.25, 0.25, nc) 0.3];
            
            lstr = '$\frac{w''}{2\epsilon\tau_o/(\rho fR)}$';
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
%     plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
grid on
set(gca, 'clim', clim);
set(gca, 'XTick', -2:1:2, 'YTick', -2:1:2);
set(gca, 'FontSize', 12);
set(get(gca, 'ylabel'), 'FontSize', 18);
% if i==1; title('Anticyclone'); end
end
set(gcf,'Color', 'w', 'Position', [384   130   613   843])


%%
ind = find(out.y == 0);

uvels = out.UMAG(:, ind);
uvels(1:ind) = -uvels(1:ind);
uvels(ind) = 0;

om = out.OMEGA(:,ind);
zetas = out.ZETA(:,ind) - om;
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

% set(gca, 'Units', 'inches');
pos = get(gca, 'position');
% set(gca, 'Position', [1.18 pos(2) 1.375 pos(4)]);
set(gca, 'Position', [.191 pos(2) 0.212 pos(4)]);
title('Anticyclone');
set(gca, 'FontSize', 12);

% legend('$\bar{u}/(\epsilon f L)$', '$\Omega/(2\epsilon f)$', '$\zeta_s/(2\epsilon f)$', 'Location', 'EastOutside');

%%
%Cyclone
out = calcCircW(epsilon, tau, f, cr);
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic;
%%
for i=1:4
    switch i
        case 1
            var = MUBAR./(epsilon*tau./f);
            conts =  linspace(-1.2, 1.2, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_x}{\tau_o/(\rho f)}$';
            pos = [4 6];
        case 2
            var = MVBAR./(tau./f);
            conts = linspace(-1.5, -.5, nc);
            conts = linspace(-1.2, -.8, nc);

            clim = conts([1 end]);
            lstr = '$\frac{M_y}{\tau_o/(\rho f)}$';
            pos = [8 10];
        case 3 
            var = W./(tau./(f*cr));
            conts = linspace(-0.25, 0.25, nc);
            clim = conts([1 end]);
            lstr = '$\frac{w_e}{\tau_o/(\rho fR)}$';
            pos = [12 14];
        case 4
            var = (W - Wclassic)./(2*epsilon*tau./(f*cr));
            clim = conts([1 end]);
            conts = [-0.3 linspace(-0.25, 0.25, nc) 0.3];
       
            lstr = '$\frac{w''}{\tau_o/(\rho fR)}$';
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
%     plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
grid on
set(gca, 'clim', clim);
set(gca, 'XTick', -2:1:2, 'YTick', -2:1:2);
set(gca, 'FontSize', 12);
if i==1;  end
end
colormap(flipud(othercolor('RdBu11')));
set(gcf,'Color', 'w', 'Position', [384   130   613   843])

%%
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
legend('$\bar{u}/(\epsilon f R)$', '$\Omega/(2\epsilon f)$', '$\zeta_s/(2\epsilon f)$', 'Location', 'EastOutside');
%% EXPORTING
% export_fig('CircularWNew.eps', '-eps', '-opengl', '-q100', '-p01')
%%
% colormap(flipud(othercolor('RdBu11')));
% clear cmnew
% cm= colormap;
% cmnew(1:32,:) = cm(1:32,:);
% cmnew(34:65,:) = cm(33:end,:);
% cmnew(33,:) = [1 1 1];
%%
% wld = wlim(end)/2;
% wconts = linspace(-wld, wld, 20);
% figure
% contourf(X./cr, Y./cr, W-Wclassic, wconts); shading interp
% axis square
% 
% uconts = -(0:.1:.6);
% % uconts = [.5 .5];
% hold on
% [c, h] = contour(X./cr, Y./cr, UMAG, uconts, 'k');
% clabel(c,h)
% hold off
% set(gca, 'xlim', [-1 1], 'ylim', [-1 1]);
% set(gca, 'clim', wlim/2);
% grid on
% title(['\epsilon_{sh} = ', num2str(epsilon), '   \epsilon_k = ', num2str(maxom./f)]);