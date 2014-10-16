clf;
hold on;

dpdt = @(p) (1 - p/10)*p;
p0 = 1;
dt = 10.^-1;
t_end = 10;
T = 0:dt:t_end;

p = @(t) 10 ./ (1+9*exp(-t));
plot(T, p(T));

p0 = 1;

eu = w01_euler( dpdt, p0, dt, t_end );
plot(T, eu);
hold off;