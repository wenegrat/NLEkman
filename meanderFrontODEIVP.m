function output = meanderFrontODEIVP(l, omega, zeta, ubar, taus, taun, f, guess)
% Numerical Solution to Meandering Front Model as IVP

OMEGA = griddedInterpolant(l, omega);
ZETA = griddedInterpolant(l, zeta);
TAUS = griddedInterpolant(l, taus);
TAUN = griddedInterpolant(l, taun);


bvpops = bvpset('NMax', 200000);

sol = ode45(@MFode, [l(1) l(end)+1], guess);
out = deval(sol, l);

output.v = out(2,:);
output.u = out(1,:);
output.l = l;

function ddl = MFode(l, y)

    omegatemp = OMEGA(l);
    zetatemp = ZETA(l);
    taustemp = TAUS(l);
    tauntemp = TAUN(l);
    r = 1e-5;
    ddl = [1./ubar*(taustemp + (f+zetatemp).*y(2)) - r.*y(1)./ubar; %Solving for u first
            1./ubar*(tauntemp - (f+2*omegatemp)*y(1))] - r.*y(2)./ubar;

end


end