function output = SolveFrontEkman(x, y, ubar, zetabase,tau,  f)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SolveFrontEkman
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Solves the Ekman flow problem (could redo MeanderIVPMulti to use this)...
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
deltax = x(2)-x(1); %assuming uniform spacing
% Calculate k
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

    % Alongfront distance
l = cumtrapz(x, abs(vels));
    
omega = k*ubar;
zeta = zetabase + omega;
     
% Project stress on BN coordinates
 taus = dot([real(tau); imag(tau)], [real(frntvec); imag(frntvec)]);
 taun = dot([real(tau); imag(tau)], [-imag(frntvec); real(frntvec)]);
     
     guess  = [0 0];
     r = 0; % Undamped in this case
     % SOLVE ODE
     out = meanderFrontODEIVP(l, omega, zeta, ubar, taus, taun, f, guess, r);
     
     % Save All Output into struct
     output.ux = real(out.u.*frntvec + out.v.*1i.*frntvec);
     output.vy =  imag(out.u.*frntvec + out.v.*1i.*frntvec);
     output.zeta = zeta;
     output.omega = omega;
     output.l =l;
     output.u = out.u;
     output.v = out.v;
     output.tangent = frntvec;
     output.taus = taus;
     output.taun = taun;
end


