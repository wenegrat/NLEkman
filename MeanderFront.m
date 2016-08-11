%Note this leads to a non-constant along-front velocity...

U0 = 1000;

r = 100;

t = 0:.01:(2*pi);

y = r.*sin(-t)
x = r.*cos(-t) + U0*t;

u = -r.*sin(-t) + U0;
v = -r.*cos(-t);
mag = abs(u+1i*v);
subplot(2,1,1)
plot(x, y, 'x')
subplot(2,1,2)
plot(u, v, 'x')

%%
hold on
for i=1:length(t)
   plot(x(i), y(i), 'x');
   pause(.001)
end
hold off

%%
