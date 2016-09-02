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

MeanderIVPMulti
%%
figure
for si=3;1:length(avecs)
    sti = find(ytemp>=avecs(si)-100,1);
    if (sti<5500); sti=5500; end
    inds = sti:(sti+2*lp);
    subplot(length(avecs), 2, 2*si-1)
    [c, h]  = contourf((x(inds)-x(inds(1)))./xfact, ys/1e3, squeeze(Wfull(4,inds,:)).', 20); shading interp
    set(h, 'edgecolor','none')
    hold on;
    contour((x(inds)-x(inds(1)))./xfact, ys/1e3, squeeze(UBfull(4,inds,:)).', uvecs, 'LineColor', 'k');
    plot((x(inds)-x(inds(1)))./xfact, ytemp(inds)./1e3, 'k')
    hold off
    ylabel(['A = ',num2str(avecs(si)/1e3), ' km']);
    set(gca, 'clim', [-1 1]*4e-5);
    set(gca, 'ylim', [-1 1]*1.5*avecs(end)/1e3)
    colorbar;
    grid on
end
% set(gcf, 'Color', 'w', 'Position', [ 675   262   905   713]);


%%
taumag = 1i*taumag;
% MERIDIONAL WIND STRESS
MeanderIVPMulti
%%
for si=1:length(avecs)
    sti = find(ytemp>=avecs(si)-100,1);
    if (sti<5500); sti=5500; end
    inds = sti:(sti+2*lp);
    subplot(length(avecs), 2, 2*si)
    [c, h]  = contourf((x(inds)-x(inds(1)))./xfact, ys/1e3, squeeze(Wfull(4,inds,:)).', 20); shading interp
    set(h, 'edgecolor','none')
    hold on;
    contour((x(inds)-x(inds(1)))./xfact, ys/1e3, squeeze(UBfull(4,inds,:)).', uvecs, 'LineColor', 'k');
    plot((x(inds)-x(inds(1)))./xfact, ytemp(inds)./1e3, 'k')
    hold off
    ylabel(['A = ',num2str(avecs(si)/1e3), ' km']);
    set(gca, 'clim', [-1 1]*4e-5);
    set(gca, 'ylim', [-1 1]*1.5*avecs(end)/1e3)
    colorbar;
    grid on
end
set(gcf, 'Color', 'w', 'Position', [ 675   262   905   713]);