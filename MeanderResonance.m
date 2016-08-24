% SOLVE VIA IVP

xfact = 120e3;
deltax = 2e3./xfact;
x = (0:deltax:500).*xfact;
f = 1e-4;
kf = 2*pi./xfact;
kn = 2*kf;
ubarc = f./kn;
taumag =0.1/1030;
epsilon = .5;
% ubarc =1;
A = 10e3;

ubar = ubarc;
dudn = 0;
% epsk = A.*(2*pi/xfact).^2./(ubar*f);


facAmp=1;

x1=x;

%Impose curvature
k = -A.*(2*pi./xfact).^2.*sin(2*pi./xfact.*x);
epsk = max(k).*ubar./f;
disp(['Eps-Omega = ', num2str(max(k).*ubar./f,2)]);


rampwidth = 100*xfact;
xc = 200.*xfact;
facTau = 1/2.*(1+tanh( (x-xc)./rampwidth));
% facTau = 1;
tau = facTau.*taumag.*ones(size(x));
taus = tau; taun=0*tau;

%%
% SOLVE ODES
omega = ubar*k;
zeta = -dudn + omega;

l = x;

guess = [0 0]; %IVP bc u, v
out = meanderFrontODEIVP(l, omega, zeta, ubar, taus, taun, f, guess);

figure 
subplot(2,1,1)
plot(x./1000, taus);
subplot(2,1,2)
plot(x./1000, out.v);



%%
disp(['L_f/L_n = ' num2str(xfact.*(f./ubar)./(2*pi),3)]);
disp(['q = 3/2\eps_k = ', num2str(3/2.*epsk)])
    

%%
% Floquet theory:
MeanderNumResonance;
disp(['Nu = ', num2str(nu)]);
