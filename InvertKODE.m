function output = InvertKODE(k,x)
% Numerical Solution to Meandering Front Model as IVP

K = griddedInterpolant(x, k);

bvpops = bvpset('NMax', 200000);
solinit = bvpinit([x(1) x(end)], [1 0 0])
sol = ode15s(@MFode, [x(1) x(end)+1], [0 -1./sqrt(2) 0]);
% sol = bvp4c(@MFode, @MFbc, solinit);

out = deval(sol, x);

output.y = out(1,:);
output.dy = out(2,:);
output.x = x;

function ddx = MFode(x, y)

   ktemp = K(x);
%    
%     ddx = [y(2); %dy1/dx = y(2)
%             ktemp.*(1+y(2).^2).^(3/2)];

   ddx=[y(2); %y(2) = TangentLine
       ktemp.*y(3); %dTdl = k*N
       -ktemp.*y(2)];
end
function res = MFbc(ya, yb)
 res = [ya(1)-yb(1) 
        ya(2)-yb(2)
        ya(3)-yb(3)];
 
end

end