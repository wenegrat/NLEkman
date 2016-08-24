%%
% Example Ekman Problem

f = 1e-4;
nits = 200;
delt = 1200;
t = (0:delt:86400*nits);
figure
hold on
omegas = [1/2 1/3 1/4 1/5 1/6].*f;
% omegas = 1.1.*omegas;
for i=1:length(omegas)
omega = omegas(i);
eps = .25;
tc = 2*pi./f .*25;
rampwidth = 2*pi./f .*20;
facAmp = 1/2.*(1+tanh( (t-tc)./rampwidth));
% facAmp = ones(size(t));
% F = sqrt(f.*(f+f.*eps.*cos(omega*t))).*t;
F = f.*(1+facAmp.*eps.*sin(omega*t)).*t;
% t2 = t(1:floor(length(t)));
t2 = (0:delt:86400*nits*2);
t2 = t;
    %%
tc = 86400*2;
rampwidth = 2*pi./f .*1;
facAmp = 1/2.*(1+tanh( (t2-tc)./rampwidth));
tau = .1.*facAmp;
% plot(tau)

% % M = - imag(exp(1i*F.*t).*cumtrapz(t, tau.*sin(-F.*t)));
% M = (conv( tau,exp(1i*F), 'full')*delt);

%%%%%%%% Matlabs trapz-function %%%%%%%%%%%%%%
M = zeros(1,numel(t));
Ft = F./t;
Ft(1) = 0;
    Fj = cumtrapz(t, Ft);
for j = 2:numel(t)
    M(j) = trapz(t(1:j), tau(1:j).*exp(-1i*(-Fj(j)+Fj(1:j))));
end

My = -imag(M);
Mx = real(M);
% plot( t(1:length(tau))./86400, M(1:length(tau)));
plot(t.*f./(2*pi), My(1:length(F))./(tau(end)./f))
% plot(t.*f, Mx(1:length(F)));
% plot(t.*f, -tau(1:length(F))./(F./t), '--');
end
hold off
%%
legend(num2str(omegas.'/f))
title(['Ekman Transport    F = f(1+\epsilon cos(\omega t))      \epsilon = ', num2str(eps)])
xlabel('Inertial Periods');
ylabel('Transport (Normalized)');
grid on