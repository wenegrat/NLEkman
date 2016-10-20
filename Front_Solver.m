function A = Front_Solver(l, xfact,ldiff, omega, zeta, f, ubar)
%Hill_Solver solves y" + w*w*f(t)*y = 0.
%Input parameters: options.w = w; - parameter 
%                  options.f = f; - function of t
%                  options.t = t; 
%f(t) is periodic, with period equal to 1.
%f(t) is normalized, such that max(f) = 1.
%The obtained solutions y_1 and y_2 are at t=1 for initial conditions (t=0)
%y_1 = 1; y'_1 = 0; y_2 = 0; y'_2 = 1.
%The output is a matrix: [y_1  y_2 ]
%                        [y'_1 y'_2]

%**************************************************************************
%**************************************************************************
OMEGA = griddedInterpolant(l, omega, 'pchip');
domegal = NaN(length(omega),1);
% domegal(1:end-1) = (omega(2:end)-omega(1:end-1))./(l(2:end)-l(1:end-1));
domegal(1:end-1) = (omega(2:end)-omega(1:end-1))./(ldiff(2:end)-ldiff(1:end-1));
domegal(end) = domegal(1);
DOMEGAL = griddedInterpolant(l, domegal,'pchip');
ZETA = griddedInterpolant(l, zeta, 'pchip');


lspan = (0:.5:1).*xfact;
[~, Y1] = ode113(@MFode, lspan, [1 0]);
A(1,1) = Y1(3,1);
A(2,1) = Y1(3,2);


[~, Y1] = ode113(@MFode, lspan, [0 1]);
A(1,2) = Y1(3,1);
A(2,2) = Y1(3,2);


function ddl = MFode(l, y)

    omegatemp = OMEGA(l);
    domegaltemp = DOMEGAL(l);
    zetatemp = ZETA(l);
%     taustemp = TAUS(l);
%     tauntemp = TAUN(l);
%     dtauntemp = DTAUN(l);
    
%     term1 = -(f+2*omegatemp)./ubar.^2.*taustemp;
    term2 = 2*domegaltemp./(f+2*omegatemp).*y(2);
    term3 = -(f+zetatemp).*(f+2*omegatemp)./ubar.^2.*y(1);
%     term3 = -(f^2 +f.*zetatemp + f*2*omegatemp)./ubar.^2.*y(1);

%     term4 = dtauntemp./ubar;
%     term5 = -tauntemp./ubar.*2*domegaltemp./(f+2*omegatemp);
    ddl = [y(2); %y(2)=dv/dl
            term2+term3];

end
end

