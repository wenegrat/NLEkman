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
wnorm = taumag./(f*L*sqrt(2));
lnorm = xfact;
nc = 40;
figure
for si=1:length(avecs)
    sti = find(ytemp>=avecs(si)-0,1);
    if (sti<5500); sti=6500; end
    if isempty(sti); sti=length(ytemp)-(2*lp+1); end
    inds = sti:(sti+2*lp);
    subtightplot(length(avecs), 1,si, gap, margh, margw)
    [c, h]  = contourf((x(inds)-x(inds(1)))./xfact, ys/lnorm, (mask(inds,:).*squeeze(Wfull(4,inds,:))).'./wnorm, nc); shading interp
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
    
    if (si==1) % PLOT FRONT STRUCTURE
        hold on
        plot(ubart/2+0.5.*ones(size(ys)), ys./lnorm, '--k', 'LineWidth', 1.25);
        hold off
    end
end
cb = colorbar('southoutside');
set(get(cb, 'xlabel'), 'String', '$\hat{w_e}$', 'Interpreter','Latex', 'FontSize', 22)
set(gcf, 'Color', 'w', 'Position', [  665   129   692   832]);
set(cb, 'Position', [  0.1988    0.1755    0.6008    0.0081]);
% set(cb, 'Ticks', -0.75:0.25:0.75);
colormap(flipud(othercolor('RdBu11')))

% EXPORT COMMAND
% export_fig('MeanderFrontNew.eps', '-eps', '-opengl', '-q100', '-p01')
%% W' PLOT
gap = [.0175 .01]; margh = .25; margw=.25;
vclassic = -taumag./(f+ZETA+.1^2.*f); %Note including small damping term.
[~, wclassic] = gradient(vclassic, deltay);
% nc = 10;
figure
avec = 6000;
si = 1;
sti = find(ytemp>=avec(si)-0,1);
    if (sti<5500); sti=6500; end
    if isempty(sti); sti=length(ytemp)-(2*lp+1); end
    inds = sti:(sti+2*lp);
    Rmin = 1./max(ktot(inds));
Lr = L./Rmin;
 for i=1:3;
     switch i
         case 1 % W NL
             var = (mask(inds,:).*squeeze(Wfull(4,inds,:))).'./wnorm;
             tstring = '$\hat{w}_{e}$';
         case 2 % W Classic
             var = wclassic(:,inds)./wnorm;
             tstring = '$\hat{w}_{SF}$';
         case 3 % W'
             wnf = (epsilon.*taumag./(f*L).*Lr);
%              wnf = wnorm./10;
             
             var = (squeeze(Wfull(4,inds,:)).' - wclassic(:,inds))./wnf;
             tstring = '$\hat{w}''$';
     end
     subtightplot(3,1,i, gap, margh, margw);
    [c, h]  = contourf((x(inds)-x(inds(1)))./xfact, ys/lnorm, var, nc); shading interp
    set(h, 'edgecolor','none')
    hold on;
    contour((x(inds)-x(inds(1)))./xfact, ys/lnorm, squeeze(UBfull(4,inds,:)).', uvecs, 'LineColor', 'k');
    plot((x(inds)-x(inds(1)))./xfact, ytemp(inds)./lnorm, 'k')
%     contour((x(inds)-x(inds(1)))./xfact, ys/lnorm, squeeze(mask(inds,:)).', [1 1], 'LineStyle', '--', 'LineColor', 'k');
    hold off
    ylabel(['A = ',num2str(avecs(si)/1e3), ' km']);
    set(gca, 'clim', [-1 1]*1);
    set(gca, 'ylim', [-1 1]*2*avecs(end)/lnorm)
    if (i~=3)
        set(gca, 'XTickLabel', []);
    else
        xlabel('$\hat{x}$', 'Interpreter', 'Latex', 'FontSize', 20);
    end
        ylabel('$\hat{y}$', 'Interpreter', 'Latex', 'FontSize', 20);
        rmin = min(abs(1./ktot(inds)));
%         t = text(.075, .075, ['R_{min} = ',num2str(rmin./lnorm, 2)]);
%             tstring = num2str(max(ytemp(inds(1)))./lnorm,2);
%             if avecs(si)==0
%                 tstring = '0';
%             end
        t = text(.075, .125, [tstring]);
        set(t, 'BackgroundColor', 'w', 'EdgeColor', 'k', 'Interpreter', 'Latex');
%         annotation('textbox',[.1 ,'String',str,'FitBoxToText','on');
        
%     colorbar;
    grid on

 end
 cb = colorbar('southoutside');
% 
set(get(cb, 'xlabel'), 'String', '$\hat{w}$', 'Interpreter','Latex', 'FontSize', 22)
% set(gcf, 'Color', 'w', 'Position', [  665   129   692   832]);
% set(cb, 'Ticks', -0.75:0.25:0.75);
set(cb, 'Position', [  0.2509    0.134    0.500    0.0158]);

colormap(flipud(othercolor('RdBu11')))
set(gcf, 'Color', 'w', 'Position', [675   377   672   587]);

% EXPORT COMMAND
% export_fig('MeanderFrontWPrime.eps', '-eps', '-opengl', '-q100', '-p01')