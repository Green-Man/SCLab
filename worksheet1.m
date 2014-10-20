clear all;
clf;
hold on;

%Time boundaries
t0 =0;dts = [1/1 1/2 1/4 1/8];t_end = 10.^1;
dts = [1/2];

%Given population ODE and initial condition
dpdt = @(p) (1 - p/10)*p;
p0 = 1;

for dt = dts
    %Execution of numerical methods
    eu = numerical( dpdt, p0, t0, dt, t_end , 'euler');
    hu = numerical( dpdt, p0, t0, dt, t_end , 'heun');
    rk = numerical( dpdt, p0, t0, dt, t_end , 'rk');
    
    T = t0:dt:t_end;
    %Analytical solution and its plot
    p = @(t) 10 ./ (1+9*exp(-t));
    P = p(T);
    
    %Plots of numerical methods
    figure(1)
    
    subplot(3,1,1);
    
%     l = legend('show');
    plot(T, P, 'Color', [0.8 0.8 0.8]);
    plot(T, eu, 'Color', [dt^0.5 0 0]);


    title('Euler')
%     s = l.TextColor;
%     l.TextColor = 'red';%{sprintf('dt = %f', dt),};
%     legend(sprintf('dt = %f', dt), 'Location','northoutside','Orientation','horizontal')
    
    subplot(3,1,2);
    h = plot(T, P, 'Color', [0.8 0.8 0.8]);
    plot(T, hu, 'Color', [0 dt^0.5 0]);
    title('Heun')
%     legend(dts)
%     
    subplot(3,1,3);
	plot(T, P, 'Color', [0.8 0.8 0.8]);
	plot(T, rk, 'Color', [0 0 dt^0.5]);
    title('Runge-Kutta')
%     legend(dts)

    

    %Plot absolute errors
%     plot(T, eu-P, 'Color', [dt^0.5 0 0]);
%     plot(T, hu-P, 'Color', [0 dt^0.5 0]);
%     plot(T, rk-P, 'Color', [0 0 1]);
    
    %Compute approximation error
    ape = @(N) sqrt(sum((N-P).^2)*dt/5);
    ape(eu);
    
end

hold off;