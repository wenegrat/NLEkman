function output = LagrangianIVP(t, guess, f, PFIELDX, PFIELDY, TAUFIELD, X, Y)
% Numerical Solution to Lagrangian equations
PX = griddedInterpolant(X.', Y.', PFIELDX.', 'linear', 'none');
PY = griddedInterpolant(X.', Y.', PFIELDY.', 'linear', 'none');
TX = griddedInterpolant(t, real(squeeze(TAUFIELD(:,1,1)))); %Assuming spatial uniform
TY = griddedInterpolant(t, imag(squeeze(TAUFIELD(:,1,1))));
bvpops = bvpset('NMax', 200000);

sol = ode45(@MFode, [t(1) t(end)+1], guess);
out = deval(sol, t);

output.v = out(4,:);
output.y = out(3,:);
output.u = out(2,:);
output.x = out(1,:);
output.t = t;

function ddt = MFode(t, y)
% 
%     omegatemp = OMEGA(l);
%     zetatemp = ZETA(l);
%     taustemp = TAUS(l);
%     tauntemp = TAUN(l);
%     
    pxtemp = PX(y(1), y(3)); 
    txtemp = TX(t);
    pytemp = PY(y(1), y(3));
    tytemp = TY(t);
    if (~isfinite(pxtemp))
%         disp('here');
    end
    ddt = [y(2); %y(1) is x so y(2) is u
        +f*y(4) + pxtemp+txtemp; %Solving for u first
           y(4); %y(3) is y so y(4) is v
            -f*y(2) + pytemp + tytemp];

end


end