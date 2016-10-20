%%% MeanderIVPMulti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MeanderIVPMulti
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
% This script is called from MeanderPlots.m
% It solves the ODEs along each path, and then interpolates to a regular
% grid.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = (0:deltax:500).*xfact; %deltax is defined in MeanderPlots.m
deltay = 500;
ys = -3*L:deltay:3*L;
ubart = ubarc.*exp(-(ys./L).^2/2); % Cross-front velocity profile

%Pre-Allocate
Wfull = NaN(length(avecs), length(x), length(ys));
UBfull = Wfull;
yfull = NaN(length(avecs), length(x));
omegatot = NaN(length(ys), length(x));
zetatot = omegatot;

%This for loop is leftover from a previous version, now solving only once
%with increasing amplitude meanders.
for Aind = 4;1:length(avecs);
    A = avecs(Aind);
    
    Mutot = NaN(length(ys), length(x));
    Mvtot = Mutot;
    ytot = Mutot;
    
    % Primary for loop (parallel for speed) that solves along each path
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    parfor ind=1:length(ys)
    ubar = ubart(ind); % Balanced vel on path
    dudn = -ubar.*ys(ind)./(L.^2); %Note assuming Gaussian vel profile
    epsk = A.*epsilon.*L.*(2*pi/xfact).^2; 
    
    % Define the meander ramp function
    rampwidth = 60*xfact;
    xc = 300.*xfact;
    facAmp = 1/2.*(1+tanh( (x-xc)./rampwidth));

    x1=x;
    offset = ys(ind);
%     y = facAmp.*A.*sin(2*pi*x./(xfact))+offset;
    y = makeOffsetCurve(facAmp.*A.*sin(2*pi*x./(xfact)), offset, x);

    ytot(ind,:) = y; % save y positions for later use

    %Determine curvature
    if (sum(~isfinite(y)) == 0)
        dx  = gradient(x, deltax*xfact); %Note that in these functions deltax is normalized by xfact
        ddx = gradient(dx, deltax*xfact); % Should be zero for uniform x
%         dy  = gradient(y, deltax*xfact);
        dy(1:end-1) = (y(2:end)-y(1:end-1))./(deltax.*xfact);
        dy(end) = dy(1);
%         ddy = gradient(dy, deltax*xfact);
        ddy(1:end-1) = (dy(2:end)-dy(1:end-1))./(deltax.*xfact);
        ddy(end) =ddy(1);
        num   = dx .* ddy - ddx .* dy; % see: http://mathworld.wolfram.com/Curvature.html Eq. (13)
        denom = dx .* dx + dy .* dy;
        denom = sqrt(denom);
        denom = denom.* denom.* denom;
        k = num ./ denom;
        k(denom < 0) = NaN;
        
        % Useful for tracking progress of parfor loop
        disp(['Eps-Omega = ', num2str(max(k).*ubar./f,2)]);

        % Determine s and n vectors
        du = gradient(x, deltax*xfact); %Temporary velocity gradients
        dv = gradient(y, deltax*xfact); %
        dv(1:end-1) = (y(2:end)-y(1:end-1))./(deltax*xfact);
        dv(end)=dv(1);
        vels = du+1i*dv;  %Tangent Vectors at each spot.
        frntvec = vels./abs(vels); %Normalized to unit length;
        
        % Define ramp function for the wind stress
        rampwidth = 20*xfact;
        xc = 40.*xfact;
        facTau = 1/2.*(1+tanh( (x-xc)./rampwidth));
        % facTau = 1;
        tau = facTau.*taumag.*ones(size(x));
        tau = tau-tau(1); %Remove first val so that it starts at exactly zero
        taus = dot([real(tau); imag(tau)], [real(frntvec); imag(frntvec)]); 
        taun = dot([real(tau); imag(tau)], [-imag(frntvec); real(frntvec)]);


        % Parameters for ODE
        omega = ubar*k;
        zeta = -dudn + omega;
        omegatot(ind,:) = omega;
        zetatot(ind,:) = zeta;
        
        % Note that solution is defined in the alongfront coordinate (not x!)
        l = abs(cumtrapz(x1, abs(vels))); % int_x sqrt( 1+(dy/dx)^2)

        ival = [0 0]; %IVP u, v
        r= 1e-5; % Define damping parameter
        
        % SOLVE ODES HERE
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        out = meanderFrontODEIVP(l, omega, zeta, ubar, taus, taun, f, ival, r); 
        
        Mu = real(out.u.*frntvec + out.v.*1i.*frntvec); %Project BN vels on cartesian coords
        Mv =  imag(out.u.*frntvec + out.v.*1i.*frntvec);

        Mutot(ind,:) = Mu; % Save for later use
        Mvtot(ind,:) = Mv;
    end
    end
    
    %% INTERPOLATE RESULTS TO CARTESIAN GRID
    npaths = length(ys); nx = length(x);
    xi = repmat(x, [1, length(ys)]).';
    yi = reshape(ytot.', npaths*nx,1);
    [X, Y] = meshgrid(x,ys); %Make new grid
    muvec = reshape(Mutot.', npaths*nx,1);
    mvvec = reshape(Mvtot.', 1,npaths*nx);
    mask = isfinite(yi); %Not all offset positions can be solved, mask those that aren't.
    
    MUBAR = griddata(xi(mask), yi(mask), muvec(mask), X, Y);
    MVBAR = griddata(xi(mask), yi(mask), mvvec(mask), X, Y);
    [MUx, MUy] = gradient(MUBAR, deltax.*xfact);
    [MVx, MVy] = gradient(MVBAR, deltay);
    W = MUx+MVy;
    Wfull(Aind, :,:) = W.';
    
    ubfull = repmat(ubart, [nx, 1]);
    ubfull = reshape(ubfull, 1, npaths*nx);
    UBAR = griddata(xi(mask), yi(mask), ubfull(mask), X, Y);
    UBfull(Aind,:,:) = UBAR.';
    zetavec = reshape(zetatot.', npaths*nx, 1);
    omegavec = reshape(omegatot.', npaths*nx, 1);
    
    OMEGA = griddata(xi(mask), yi(mask), omegavec(mask), X, Y);
    ZETA = griddata(xi(mask), yi(mask), zetavec(mask), X, Y);
    
end