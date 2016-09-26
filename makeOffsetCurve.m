function y = makeOffsetCurve(yorig, offset, x)
%% Generate an offset curve of constant offset, return at original x positions.
xdiff = x(2)-x(1);
yprime = gradient(yorig, xdiff);
xprime = 1; %Assumed input 

xs = x + offset.*yprime./(sqrt(1+yprime.^2));
ys = yorig - offset./(sqrt(1+yprime.^2));

y = interp1(xs, ys, x, 'pchip');
end