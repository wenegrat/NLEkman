function output = meanderFrontODEs(l, omega, zeta, ubar, taus, taun, f, guessVec)
% Numerical Solution to Meandering Front Model
% ubar*dv/dl + (1+2*omega)u = taun
% ubar*du/dl -(1+zeta)v = taus;
% y1 = u
% y2 = v;

% y1' = 1/ubar*[taus + (1+zeta)y2]
% y2' = 1/ubar*[taun - (1+2*omega)y1]

% UINIT = griddedInterpolant(l, guess(1,:));
% VINIT = griddedInterpolant(l, guess(2,:));
dudn = zeta - omega;
solinit.x=l;
solinit.y=guessVec;

for i=1
    omegat = omega./i;
    zeta = dudn + omegat;
OMEGA = griddedInterpolant(l, omegat);
ZETA = griddedInterpolant(l, zeta);
TAUS = griddedInterpolant(l, taus);
TAUN = griddedInterpolant(l, taun);
bvpops = bvpset('NMax', 20000, 'AbsTol', 0.001);
% solinit = bvpinit(l,@MFinit);
% solinit = bvpinit(l,@MFinit);
% solinit = bvpinit(l, guessVec);

sol = bvp4c(@MFode, @MFbc,solinit, bvpops);
out = deval(sol, l);
solinit.y = out;
end


output.u = out(1,:);
output.v = out(2,:);
% output.uz = real(out(1,:)./(1023*av));
% output.vz = imag(out(1,:)./(1023*av));
output.l = l;

%====================
function init = MFinit(l)
 init = [TAUN(l)./(f+2*OMEGA(l)) 
        -TAUS(l)./(f+ZETA(l))].';
%  init = [0; 0];
end
%===================

function ddl = MFode(l, y)

    omegatemp = OMEGA(l);
    zetatemp = ZETA(l);
    taustemp = TAUS(l);
    tauntemp = TAUN(l);
    r = 0;1e-5;
    ddl = [1./ubar*(taustemp + (f+zetatemp).*y(2)) - r.*y(1)./ubar; %Solving for u first
            1./ubar*(tauntemp - (f+2*omegatemp)*y(1))] - r.*y(2)./ubar;

end

%====================
function res = MFbc(ya, yb)
 res = [ya(1)-yb(1) 
        ya(2)-yb(2)];
 
end

end