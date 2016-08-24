function out = pathToGrid(input, curves,x, y)
[nc, nx] = size(curves);

out = zeros(nx, length(y));
for i=1:nc
    for j=1:length(x)
        
    out(j,:) = interp1(curves(:,j), input(:,j), y);
    end
end
end