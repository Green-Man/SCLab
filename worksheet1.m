clear all;
clf;
hold on;

%Time boundaries
t0 =0;
dts = [1/1 1/2 1/4 1/8];
t_end = 10.^1;


%Given population ODE and initial condition
dpdt = @(p) (1 - p/10)*p;
p0 = 1;

methods = {'euler', 'heun', 'rk'};
i = 1;
for method = methods
    figure(1)
    hold on
    subplot(size(methods,2),1,i);
    legendStrings = {'Analytic'};
    for dt = dts

        T = t0:dt:t_end;
        %Analytical solution and its plot
        p = @(t) 10 ./ (1+9*exp(-t));
        P = p(T);
        plot(T, P, 'Color', [0.8 0.8 0.8]);
        
        %Execution of numerical methods
        ns = numerical( dpdt, p0, t0, dt, t_end , method);
        
        %Plots of numerical methods
        gColor = [0,0,0]; gColor(i) = dt^0.5;
        hold on;
        plot(T, ns, 'Color', gColor);
        hold on;
        
        
    %     s = l.TextColor;
    %     l.TextColor = 'red';%{sprintf('dt = %f', dt),};
    %     legend(sprintf('dt = %f', dt), 'Location','northoutside','Orientation','horizontal')

        %Plot absolute errors
%         plot(T, ns-P, 'Color', [dt^0.5 0 0]);

        %Compute approximation error
        ape = @(N) sqrt(sum((N-P).^2)*dt/5);
        legendStrings = [legendStrings, {sprintf('dt=%2.3f ae=%2.6f', dt, ape(ns))}];
    end
    
    legend(legendStrings)
    title(method);
    i=i+1;
end

hold off;