function [omega, phi, ub, vb] =calculateOmega(ETA,U, V, dx, dy)
g = 9.81;

f = 1e-4;

[dEtady, dEtadx, ~, ~] = gradient(ETA, dx);
% Find balanced vel field
ub = -g./f.*dEtady;
vb =  g./f.*dEtadx;

ub = U(:,:,end,:);
vb = V(:,:,end,:);

% calculate phi angle
phi = angle(ub+1i.* vb);
% Calculate Omega
[dPhidy, dPhidx, ~, ~] = gradient(phi, dx);

[nx, ny, nz, nt] = size(dPhidy);
dvec = reshape(dPhidx, nx*ny*nz*nt,1);
st = nanstd(dvec);
dPhidx(abs(dPhidx)>2*st) = NaN;
[xx, yy, tt] = ndgrid(1:nx, 1:ny, 1:nt);
mask = isfinite(squeeze(dPhidx));
dp = squeeze(dPhidx);
dPhidx(:,:,1,:) = interpn(xx, yy, tt, dp, xx, yy, tt);
omega = ub.*dPhidx + vb.*dPhidy;


tvec = (ub + 1i*vb)./(abs(ub+1i*vb));
nvec = -imag(tvec)+1i*real(tvec);

[uy, ux] = gradient(ub, dx);
[vy, vx] = gradient(vb, dx);

mag = abs(ub+1i*vb);
[magy, magx] = gradient(mag, dx);

zeta = vx - uy;
zetas = - (real(nvec).*magx + imag(nvec).*magy);
omega = zeta - zetas;

% ts = 10800;
% 
% du  = gradient(ub*ts, dx);
% [~, ddx] = gradient(du, dx);
% dv  = gradient(vb*ts, dx);
% [ddy, ~] = gradient(dv, dx);
% 
% num   = du .* ddy - ddx .* dv;
% denom = du .* du + dv .* dv;
% denom = sqrt(denom);
% denom = denom.* denom.* denom;
% ktot = num ./ denom;
% ktot(denom < 0) = NaN;
% 
% omega = abs(ub+1i*vb).*ktot;
end