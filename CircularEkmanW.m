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

xl = [-1.5 1.5];

gap = [.03 .01]; margh = .1; margw = .1;

for i=1:4
    switch i
        case 1
            var = MUBAR./(epsilon.*tau./f);
            conts =  linspace(-2.5, 2.5, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_x}{\epsilon\tau_o/(\rho f)}$';
        case 2
            var = MVBAR./(tau./f);
            conts = linspace(-1.5, -.5, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_y}{\tau_o/(\rho f)}$';
        case 3 
            var = W./(tau./(f*cr));
            conts = linspace(-1.5, 1.5, nc);
            clim = conts([1 end]);
            lstr = '$\frac{w_e}{\tau_o/(\rho fR)}$';
        case 4
            var = (W - Wclassic)./(tau./(f*cr));
            conts = conts./5;
            clim = conts([1 end]);
            lstr = '$\frac{w''}{\epsilon\tau_o/(\rho fR)}$';
%             conts = wcontsclass;
    end
subtightplot(4,2,2*i-1, gap, margh, margw);
[c, h] = contourf(X./cr, Y./cr, var, conts); shading interp
% set(h, 'edgecolor', 'none');
axis square
% colorbar;
set(gca, 'xlim', xl, 'ylim', xl);
% grid on
ylabel(lstr, 'Rotation', 0, 'Interpreter', 'Latex');
set(get(gca, 'ylabel'), 'Position', [-2.6468   -0.2819         0]);
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
    plot(cos(0:.1:2*pi), sin(0:.1:2*pi), '--k'); 
    plot(0.65.*cos(0:.1:2*pi), 0.65.*sin(0:.1:2*pi), '--k'); 
%     plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
grid on
set(gca, 'clim', clim);
set(gca, 'FontSize', 12);
set(get(gca, 'ylabel'), 'FontSize', 18);
if i==1; title('Anticyclone'); end
end


%%
%Cyclone
out = calcCircW(epsilon, tau, f, cr);
X = out.x; Y = out.y; MUBAR = out.MU; MVBAR = out.MV; W = out.W; Wclassic = out.WClassic;
%%
for i=1:4
    switch i
        case 1
            var = MUBAR./(epsilon.*tau./f);
            conts =  linspace(-2.5, 2.5, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_x}{\epsilon\tau_o/(\rho f)}$';
        case 2
            var = MVBAR./(tau./f);
            conts = linspace(-1.5, -.5, nc);
            clim = conts([1 end]);
            lstr = '$\frac{M_y}{\tau_o/(\rho f)}$';
        case 3 
            var = W./(tau./(f*cr));
            conts = linspace(-1.5, 1.5, nc);
            clim = conts([1 end]);
            lstr = '$\frac{w_e}{\tau_o/(\rho fR)}$';
        case 4
            var = (W - Wclassic)./(tau./(f*cr));
            conts = conts./5;
            clim = conts([1 end]);
            lstr = '$\frac{w''}{\epsilon\tau_o/(\rho fR)}$';
%             conts = wcontsclass;
    end
subtightplot(4,2,2*i, gap, margh, margw);
[c, h] = contourf(X./cr, Y./cr, var, conts); shading interp
% set(h, 'edgecolor', 'none');
axis square
cb = colorbar;
set(gca, 'xlim', xl, 'ylim', xl);
% grid on
% ylabel(lstr, 'Rotation', 0, 'Interpreter', 'Latex');
% hold on; [c, h] = contour(X./cr, Y./cr, UMAG, uconts, '-k'); hold off
hold on; 
    plot(cos(0:.1:2*pi), sin(0:.1:2*pi), '--k'); 
    plot(0.65.*cos(0:.1:2*pi), 0.65.*sin(0:.1:2*pi), '--k'); 
%     plot(0.25.*cos(0:.1:2*pi), 0.25*sin(0:.1:2*pi), 'k'); 
hold off
grid on
set(gca, 'clim', clim);
set(gca, 'FontSize', 12);
if i==1; title('Cyclone'); end
end
colormap(flipud(othercolor('RdBu11')));
set(gcf,'Color', 'w', 'Position', [384   130   613   843])

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