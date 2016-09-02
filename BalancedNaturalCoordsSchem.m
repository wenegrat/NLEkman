deltax  = .1;
x = 0:deltax:2*pi+.1;

y = 1*sin(x);
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
    vels = dx+1i*dy;  %Tangent Vectors at each spot.
    frntvec = vels./abs(vels); %Normalized;
    normal = 1i*frntvec;
    
    pointind = 15;
    point = x(pointind)+1i*y(pointind);
    R = 1./k(pointind);
    normcirc = normal(pointind);
%     if (R<0) normcirc = -normcirc; end
    circcent = point + R.*normcirc;
    
    vecind = 45;
    
    figure
plot(x,y, 'LineWidth', 2, 'Color', 'k');
hold on
plot(x(pointind), y(pointind), 'kx');
plot(R*cos(x)+real(circcent), R*sin(x)+imag(circcent),'k');
h = line([x(pointind) real(circcent)], [y(pointind) imag(circcent)]);
set(h, 'Color', 'k');
% annotation('arrow',[x(end-1) x(end)], [y(end-1) y(end)])
annotation('arrow',[0.450877192982456 0.507017543859649],...
    [0.600895734597156 0.533175355450237]);
annotation('arrow',[0.856140350877193 0.889473684210526],...
    [0.446867298578199 0.502369668246446]);
h = line([x(vecind) x(vecind)+real(frntvec(vecind))], [y(vecind) y(vecind)+imag(frntvec(vecind))]);
set(h, 'Color', 'r');
h = line([x(vecind) x(vecind)-imag(frntvec(vecind))], [y(vecind) y(vecind)+real(frntvec(vecind))]);
set(h, 'Color', 'r');
hold off
grid on;
axis equal
set(gca, 'xlim', [0 2*pi], 'ylim', [-1.5 1.5]);
set(gca, 'XTickLabel', [], 'YTickLabel',[]);