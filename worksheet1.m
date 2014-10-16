clf;
hold on;

%Time boundaries
t0 =0;dt = 10.^-0.3;t_end = 10.^1;
T = t0:dt:t_end;
iter_to_plot = round(t_end/dt);

%Analytical solution and its plot
p = @(t) 10 ./ (1+9*exp(-t));
plot(T(1:iter_to_plot), p(T(1:iter_to_plot)), 'Color', [0.8 0.8 0.8]);

%Given population ODE and initial condition
dpdt = @(p) (1 - p/10)*p;
p0 = 1;

%Execution of numerical methods
eu = w01_euler( dpdt, p0, t0, dt, t_end );
hu = w01_heun( dpdt, p0, t0, dt, t_end );
rk1 = w01_rk1( dpdt, p0, t0, dt, t_end );

%Plots of numerical methods
plot(T(1:iter_to_plot), eu(1:iter_to_plot), 'Color', [1 0 0]);
plot(T(1:iter_to_plot), hu(1:iter_to_plot), 'Color', [0 1 0]);
plot(T(1:iter_to_plot), rk1(1:iter_to_plot), 'Color', [0 0 1]);
hold off;