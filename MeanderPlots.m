%% 
% PLOTS FOR MEANDERIVP
f=1e-4;
taumag = .1./1035;
avecs = [0 3e3  6e3 12.25e3]; %Note this is no longer relevant, doing only 1 integration.
xfact = 120e3;
deltax = (2000./xfact);
% ZONAL WIND STRESS
%%

lp = xfact./(deltax.*xfact);
epsilon = .5;
ubarc = .5;
L = ubarc./(f*epsilon);
uvecs = [.3679 .3679].*ubarc;

%CALCULATE SOLUTION
MeanderIVPMulti;

%%
rampwidth = 60*xfact;
xc = 300.*xfact;
facAmp = 1/2.*(1+tanh( (x-xc)./rampwidth));
ytemp = avecs(end).*facAmp.*sin(2*pi./xfact.*x);
    dx  = gradient(x, deltax*xfact);
    ddx = gradient(dx, deltax*xfact);
    dy  = gradient(ytemp, deltax*xfact);
    dy(1:end-1) = (ytemp(2:end)-ytemp(1:end-1))./(deltax.*xfact);
    dy(end) = dy(1);
    ddy = gradient(dy, deltax*xfact);
    ddy(1:end-1) = (dy(2:end)-dy(1:end-1))./(deltax.*xfact);
    ddy(end) =ddy(1);
    num   = dx .* ddy - ddx .* dy;
    denom = dx .* dx + dy .* dy;
    denom = sqrt(denom);
    denom = denom.* denom.* denom;
    ktot = num ./ denom;
    ktot(denom < 0) = NaN;
    
%%
mask = ones(size(Y));
[ny nx] = size(Y);
for i=1:nx
    maxoffset = find(isfinite(ytot(:,i)), 1, 'last');
    minoffset = find(isfinite(ytot(:,i)), 1, 'first');
    for j=1:ny
       
       mask(j,i) = (Y(j,i)>=ytot(minoffset,i))& ( Y(j,i)<= ytot(maxoffset,i));
    end
end
mask = mask.';
%%
avecs = [0 3e3 6e3 12e3];
gap = [.02 .02]; margh=.25; margw=.2;
wnorm = taumag./(1+epsilon).^2.*epsilon; % XX-Check this.
wnorm = taumag./(f*L);
lnorm = xfact;
figure
for si=1:length(avecs)
    sti = find(ytemp>=avecs(si)-0,1);
    if (sti<5500); sti=6500; end
    if isempty(sti); sti=length(ytemp)-(2*lp+1); end
    inds = sti:(sti+2*lp);
    subtightplot(length(avecs), 1,si, gap, margh, margw)
    [c, h]  = contourf((x(inds)-x(inds(1)))./xfact, ys/lnorm, (mask(inds,:).*squeeze(Wfull(4,inds,:))).'./wnorm, 20); shading interp
    set(h, 'edgecolor','none')
    hold on;
    contour((x(inds)-x(inds(1)))./xfact, ys/lnorm, squeeze(UBfull(4,inds,:)).', uvecs, 'LineColor', 'k');
    plot((x(inds)-x(inds(1)))./xfact, ytemp(inds)./lnorm, 'k')
%     contour((x(inds)-x(inds(1)))./xfact, ys/lnorm, squeeze(mask(inds,:)).', [1 1], 'LineStyle', '--', 'LineColor', 'k');
    hold off
    ylabel(['A = ',num2str(avecs(si)/1e3), ' km']);
    set(gca, 'clim', [-1 1]*1);
    set(gca, 'ylim', [-1 1]*2.25*avecs(end)/lnorm)
    if (si~=length(avecs))
        set(gca, 'XTickLabel', []);
    else
        xlabel('$\hat{x}$', 'Interpreter', 'Latex', 'FontSize', 20);
    end
        ylabel('$\hat{y}$', 'Interpreter', 'Latex', 'FontSize', 20);
        rmin = min(abs(1./ktot(inds)));
%         t = text(.075, .075, ['R_{min} = ',num2str(rmin./lnorm, 2)]);
            tstring = num2str(max(ytemp(inds(1)))./lnorm,2);
            if avecs(si)==0
                tstring = '0';
            end
        t = text(.075, .15, ['$\hat{A} =$ ',tstring]);
        set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
%         annotation('textbox',[.1 ,'String',str,'FitBoxToText','on');
        
%     colorbar;
    grid on
end
cb = colorbar('southoutside');
set(get(cb, 'xlabel'), 'String', '$\hat{w_e}$', 'Interpreter','Latex', 'FontSize', 22)
set(gcf, 'Color', 'w', 'Position', [  665   129   692   832]);
set(cb, 'Position', [  0.1988    0.1755    0.6008    0.0081]);
colormap(flipud(othercolor('RdBu11')))

% EXPORT COMMAND
% export_fig('MeanderFrontNew.eps', '-eps', '-opengl', '-q100', '-p01')
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