function output = meanderFrontODESep(l, omega, zeta, ubar, taus, taun, f)
% Numerical Solution to Meandering Front Model
% ubar*dv/dl + (1+2*omega)u = taun
% ubar*du/dl -(1+zeta)v = taus;
% y1 = u
% y2 = v;

% y1' = 1/ubar*[taus + (1+zeta)y2]
% y2' = 1/ubar*[taun - (1+2*omega)y1]
lover = l(end);
% UINIT = griddedInterpolant(l, guess(1,:));
% VINIT = griddedInterpolant(l, guess(2,:));
OMEGA = griddedInterpolant(l, omega);
domegal = NaN(length(omega),1);
domegal(1:end-1) = (omega(2:end)-omega(1:end-1))./(l(2:end)-l(1:end-1));
domegal(end) = domegal(1);
% domegal = omega(gradient(omega, l);
DOMEGAL = griddedInterpolant(l, domegal);
ZETA = griddedInterpolant(l, zeta);
TAUS = griddedInterpolant(l, taus);
TAUN = griddedInterpolant(l, taun);
dtaun = gradient(taun,l);
dtaun(1:end-1) = (taun(2:end)-taun(1:end-1))./(l(2:end)-l(1:end-1));
dtaun(end) = dtaun(1);
DTAUN = griddedInterpolant(l,dtaun);
bvpops = bvpset('NMax', 200000);
% solinit = bvpinit(l,@MFinit);
solinit = bvpinit(l,@MFinit);

% sol = bvp4c(@MFode, @MFbc,solinit);
 sol = ode45(@MFode, [l(1) l(end)+1], [0 0]);
out = deval(sol, l);

output.v = out(1,:);
u = 1./(f+2*omega).*(taun - ubar.*out(2,:));
output.u = u;
% output.uz = real(out(1,:)./(1023*av));
% output.vz = imag(out(1,:)./(1023*av));
output.l = l;

%====================
function init = MFinit(l)
 init = [TAUN(l)./(f+ZETA(l)) 
        -TAUS(l)./(f+2*OMEGA(l))].';
 init = [sin(2*pi*l./(lover)); cos(2*pi*l./lover)];
 init = [-.6 0];
end
%===================

function ddl = MFode(l, y)

    omegatemp = OMEGA(l);
    domegaltemp = DOMEGAL(l);
    zetatemp = ZETA(l);
    taustemp = TAUS(l);
    tauntemp = TAUN(l);
    dtauntemp = DTAUN(l);
    
    term1 = -(f+2*omegatemp)./ubar.^2.*taustemp;
    term2 = 2*domegaltemp./(f+2*omegatemp).*y(2);
%      term2 = 0;
    term3 = -(f+zetatemp).*(f+2*omegatemp)./ubar.^2.*y(1);
    term4 = dtauntemp./ubar;
    term5 = -tauntemp./ubar.*2*domegaltemp./(f+2*omegatemp);
%     term5=0;
    ddl = [y(2); %y(2)=dv/dl
            term1+term2+term3+term4+term5];

end

%====================
function res = MFbc(ya, yb)
 res = [ya(1)-yb(1) 
        ya(2)-yb(2)];
 
end

end