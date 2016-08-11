function output = meanderFrontODEs(l, omega, zeta, ubar, taus, taun, f)
% Numerical Solution to Meandering Front Model
% ubar*dv/dl + (1+2*omega)u = taun
% ubar*du/dl -(1+zeta)v = taus;
% y1 = u
% y2 = v;

% y1' = 1/ubar*[taus + (1+zeta)y2]
% y2' = 1/ubar*[taun - (1+2*omega)y1]


OMEGA = griddedInterpolant(l, omega);
ZETA = griddedInterpolant(l, zeta);
TAUS = griddedInterpolant(l, taus);
TAUN = griddedInterpolant(l, taun);
bvpops = bvpset('NMax', 200000);
% solinit = bvpinit(l,@MFinit);
solinit = bvpinit([l(1) l(end)],[1 1]);

sol = bvp4c(@MFode, @MFbc,solinit, bvpops);

out = deval(sol, l);

output.u = out(1,:);
output.v = out(2,:);
% output.uz = real(out(1,:)./(1023*av));
% output.vz = imag(out(1,:)./(1023*av));
output.l = l;

%====================
function init = MFinit(l)
 init = [taun./(f+zeta) 
        -taus./(f+2*omega)].';
end
%===================

function ddl = MFode(l, y)

    omegatemp = OMEGA(l);
    zetatemp = ZETA(l);
    taustemp = TAUS(l);
    tauntemp = TAUN(l);
    
    ddl = [1./ubar*(taustemp + (f+zetatemp).*y(2));
            1./ubar*(tauntemp - (f+2*omegatemp)*y(1))];

end

%====================
function res = MFbc(ya, yb)
 res = [ya(1)-yb(1) 
        ya(2)-yb(2)];
 
end

end