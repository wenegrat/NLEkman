kx = linspace(1e-8, 2e-4, 1000);
kxn = kx./1e-4; %kx/f

h = linspace(0, 1, 1000); % h is 3*omega./f;

% resonance when kx + eps is 2f/n for n=1, 2, 3, ...
% Landau page 84, delta eps = n^(2n-3)h^n f / (2^(3(n-1) *( (n-1)! )^2

[K, H] = meshgrid(kxn, h);
stability = zeros(size(K));
for n = 1:6;
   
    wn = 2/n;
    deltaE = n^(2*n-3).*H.^n .* wn / (2^(3*(n-1)).*( factorial(n-1).^2));
    xl = wn - deltaE;
    xr = wn + deltaE;
    criteria =  (xl <= K) & (K <= xr);
    stability(criteria) = 1;
    
end

% for j=1:length(kxn)
%    n = kxn./ 
% end

pcolor(kxn, h, stability); shading interp
colormap(flipud(gray))

%%
% Do the finite difference version (Following Kutz section 7.8).

deltaL = 1e3;
xl = 1;
N = 1000;
x = linspace(1, 5*xl, N+2);
deltaL = x(2)-x(1);
kn = 2*pi./(xl);
eps = .8;

% N = length(x);
B = zeros(N,N);
for j=1:N
    B(j,j) = -2; %Diagonal Element
end
for j=1:N-1;
    B(j,j+1)=1; %Off diagonal elements
    B(j+1,j)=1;
end
B1=B./deltaL.^2; % Dirichlet Matrix 

B3 = B;
B3(N,1) =1; B3(1,N)=-1;
B3 = B3./deltaL.^2; %Periodic matrix

% Define ODE as:
% B1*y + f^2(1+eps*cos(kx*l))y = 0
potfunc = kn.^2.*(1+eps.*cos(2*kn*x));
P=zeros(N,N);
for j=1:N
    P(j,j) = potfunc(j+1);
end

linL = B3 + P;
d = sort(eig(linL));

% d = diag(D);
plot(d);