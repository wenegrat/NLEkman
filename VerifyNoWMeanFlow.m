lambda = 120e3;
k = 2*pi./lambda;
deltax = floor(lambda./500);
x = 0:deltax:lambda;
A = 12e3;
offsets = -3*A:deltax:3*A;
ubar = 1;
y = NaN(length(offsets), length(x));
xs = y;
Mu = y;
Mv = y;
thetas = y;
for i=1:length(offsets)
    off = offsets(i);
    xs(i,:) = x + off.*A*k*cos(k*x)./(sqrt(1+k.^2.*A.^2.*cos(k*x).^2));
    y(i,:) = A*sin(k*x) - off./sqrt(1+k^2.*A^2.*cos(k*x).^2);
    xs(i,:) = x;
    y(i,:) = makeOffsetCurve(A*sin(k*x), off, x);
%     xs(i,:) = x;
%     y(i,:) = A*sin(k*x) + off;
%     u(i,:) = ubar.*ones(size(y));
    
    % Determine s and n vectors
    du = gradient(xs(i,:), deltax); %Temporary velocity gradients
    dv = gradient(y(i,:), deltax); %
    dv(1:end-1) = (y(i,2:end)-y(i,1:end-1))./(deltax);
    dv(end)=dv(1);
    vels = du+1i*dv;  %Tangent Vectors at each spot.
    frntvec = vels./abs(vels); %Normalized;
    velst(i,:) = vels;
    Mu(i,:) = real(ubar*frntvec + 0*1i.*frntvec);
    Mv(i,:) =  imag(ubar*frntvec + 0*1i.*frntvec);
    thetas(i,:) = angle(frntvec);
    
     %Determine curvature
    dx  = gradient(xs(i,:), deltax);
    ddx = gradient(dx, deltax);
%     dy  = gradient(y, deltax);
    dy = NaN(1,length(y));
    dy(1:end-1) = (y(i,2:end)-y(i,1:end-1))./(deltax);
    dy(end) = dy(1);
    ddy = gradient(dy, deltax);
    ddy(1:end-1) = (dy(2:end)-dy(1:end-1))./(deltax);
    ddy(end) =ddy(1);
    num   = dx .* ddy - ddx .* dy;
    denom = dx .* dx + dy .* dy;
    denom = sqrt(denom);
    denom = denom.* denom.* denom;
    kte = num ./ denom;
    kte(denom < 0) = NaN;
    kt(i,:) = kte;
end
   npaths = length(offsets); nx = length(x);
    xi = repmat(x, [1, length(offsets)]).';
    xi = reshape(xs.', npaths*nx, 1);
    yi = reshape(y.', npaths*nx,1);
    [X, Y] = meshgrid(x,offsets);
    muvec = reshape(Mu.', npaths*nx,1);
    mvvec = reshape(Mv.', npaths*nx,1);
    MUBAR = griddata(xi, yi, muvec, X, Y);
    MVBAR = griddata(xi, yi, mvvec, X, Y);
    [MUx, MUy] = gradient(MUBAR, deltax);
    [MVx, MVy] = gradient(MVBAR, deltax);
    
    W = MUx+MVy;
    thetavec = reshape(thetas.', npaths.*nx,1);
    THETA = griddata(xi,yi, thetavec, X, Y);
    THETA = angle(MUBAR + 1i*MVBAR);
    UTEST = ubar.*cos(THETA);
    VTEST = ubar.*sin(THETA);
    [MUxt, ~] = gradient(UTEST, deltax);
    [~, MVyt] = gradient(VTEST, deltax);
    Wt = MUxt + MVyt;
    %%
%     contourf(x, offsets, W, 100); shading interp
    pcolor(x, offsets, W); shading interp
%     contourf(xs, y, thetas, 100); 
    colorbar;
    hold on
    inds = 145:5:155;
    plot(xs(inds,:).', y(inds,:).','k');
    hold off
        axis equal
    set(gca, 'ylim', [-1 1].*A*1.25);
    set(gca, 'clim', [-1 1].*1e-6);
%    set(gca, 'clim', [-1 1].*.5);