%% 
% PLOTS FOR MEANDERIVP
taumag = .1./1035;
% ZONAL WIND STRESS
%%

lp = xfact./(deltax.*xfact);
uvecs = [.3679 .3679].*ubarc;

rampwidth = 60*xfact;
xc = 280.*xfact;
facAmp = 1/2.*(1+tanh( (x-xc)./rampwidth));
ytemp = avecs(end).*facAmp.*sin(2*pi./xfact.*x);

% MeanderIVPMulti
%%
gap = [.02 .02]; margh=.25; margw=.1;
wnorm = taumag./(1+epsilon).^2.*epsilon; % XX-Check this.
lnorm = xfact;
figure
for si=1:length(avecs)
    sti = find(ytemp>=avecs(si)-300,1);
    if (sti<5500); sti=5000; end
    if isempty(sti); sti=length(ytemp)-(2*lp+1); end
    inds = sti:(sti+2*lp);
    subtightplot(length(avecs), 1,si, gap, margh, margw)
    [c, h]  = contourf((x(inds)-x(inds(1)))./xfact, ys/lnorm, squeeze(Wfull(4,inds,:)).'./wnorm, 20); shading interp
    set(h, 'edgecolor','none')
    hold on;
    contour((x(inds)-x(inds(1)))./xfact, ys/lnorm, squeeze(UBfull(4,inds,:)).', uvecs, 'LineColor', 'k');
    plot((x(inds)-x(inds(1)))./xfact, ytemp(inds)./lnorm, 'k')
    hold off
    ylabel(['A = ',num2str(avecs(si)/1e3), ' km']);
    set(gca, 'clim', [-1 1]*4);
    set(gca, 'ylim', [-1 1]*1.5*avecs(end)/lnorm)
    if (si~=length(avecs))
        set(gca, 'XTickLabel', []);
    else
        xlabel('$\hat{x}$', 'Interpreter', 'Latex', 'FontSize', 20);
    end
        ylabel('$\hat{y}$', 'Interpreter', 'Latex', 'FontSize', 20);
        rmin = min(abs(1./ktot(inds)));
        t = text(.075, .075, ['R_{min} = ',num2str(rmin./lnorm, 2)]);
        set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k');
%         annotation('textbox',[.1 ,'String',str,'FitBoxToText','on');
        
%     colorbar;
    grid on
end
cb = colorbar('southoutside');
set(get(cb, 'xlabel'), 'String', '$\hat{w_e}$', 'Interpreter','Latex', 'FontSize', 22)
set(gcf, 'Color', 'w', 'Position', [  670   147   722   825]);
set(cb, 'Position', [0.0997    0.1685    0.8006    0.010]);

%%
% taumag = 1i*taumag;
% % MERIDIONAL WIND STRESS
% MeanderIVPMulti
% %%
% for si=1:length(avecs)
%     sti = find(ytemp>=avecs(si)-100,1);
%     if (sti<5500); sti=5500; end
%     inds = sti:(sti+2*lp);
%     subplot(length(avecs), 2, 2*si)
%     [c, h]  = contourf((x(inds)-x(inds(1)))./xfact, ys/1e3, squeeze(Wfull(4,inds,:)).', 20); shading interp
%     set(h, 'edgecolor','none')
%     hold on;
%     contour((x(inds)-x(inds(1)))./xfact, ys/1e3, squeeze(UBfull(4,inds,:)).', uvecs, 'LineColor', 'k');
%     plot((x(inds)-x(inds(1)))./xfact, ytemp(inds)./1e3, 'k')
%     hold off
%     ylabel(['A = ',num2str(avecs(si)/1e3), ' km']);
%     set(gca, 'clim', [-1 1]*4e-5);
%     set(gca, 'ylim', [-1 1]*1.5*avecs(end)/1e3)
%     colorbar;
%     grid on
% end
% set(gcf, 'Color', 'w', 'Position', [ 675   262   905   713]);