clear all;clf;hold on;
disp('============================================')
figure(1)
hold on
format shortG

%Time boundaries
t_start = 0;
dts = ones(1,6);
for dti = 1:length(dts)
    dts(dti) = 1./(2.^(dti-1));
end
dts = fliplr(dts);
% dts = dts(3);
t_end = 5*10^0;
t_end_plot = 1.5;
newtonE = 1e-4;

%Given ODE and initial condition
dpdt = @(p) 7*(1 - p/10)*p;
p0 = 20;

%Analytical solution
p = @(t) 200 ./ (20-10*exp(-7*t));
exactT = t_start:dts(1):t_end;
P = p(exactT);

methods = {'Explicit Euler', 'Heun', 'Runge-Kutta'};
methods = {'Explicit Euler', 'Heun', 'Implicit Euler'};
methodIdx = 1;
exactError = zeros(size(methods,2), size(dts,2));
approxError = zeros(size(methods,2), size(dts,2));
times = zeros(size(methods,2), size(dts,2));

for method = methods
    %Plot each method on the separate subplot
    subplot(1,size(methods,2),methodIdx);
    legendStrings = {'Analytic'};
    plot(exactT(1:round(t_end_plot/dts(1))), P(1:round(t_end_plot/dts(1))), 'Color', [0.3 0.3 0.3]);
    
    nsBest = [];	%the most precise solution
    timestepIdx=1;
    for dt = dts
        T = t_start:dt:t_end;
        
        %Execution of numerical methods
        [ns times(methodIdx,timestepIdx)] = numerical( dpdt, p0, t_start, dt, t_end , method, newtonE);
        
        %Save the most precise solution
        if dt == dts(1) nsBest = ns; end
        
        %Plot of the numerical methods
        hold on; ylim([0 20]);
        gColor = [0,0,0]; gColor(methodIdx) = min([dt^0.5 1]);
        plot(T(1:round(t_end_plot/dt)), ns(1:round(t_end_plot/dt)), 'x', 'Color', gColor);

        %Compute the exact error
        %interpP = interp1(exactT, P, t_start:dt:t_end, 'linear');
        interpP = P(1:dt/dts(1):length(P));
        exactErrorF = @(N) sqrt(sum((N-interpP).^2)*dt/5);
        exactError(methodIdx,timestepIdx) = exactErrorF(ns);
        
        %Compute approximation error
        %interpNsBest = interp1(exactT, nsBest, t_start:dt:t_end, 'linear');
        interpNsBest = nsBest(1:dt/dts(1):length(nsBest));
        approximationErrorF = @(N) sqrt(sum((N-interpNsBest).^2)*dt/5);
        if dt ~= dts(1)
            approxError(methodIdx,timestepIdx) = approximationErrorF(ns);
        end
        
        %Format a legend
        errorString = sprintf('exactError=%2.6f', exactError(methodIdx,timestepIdx));
        dtString = sprintf('dt=%2.3f', dt);
        nextLegendString = sprintf('%s', dtString);
        legendStrings = [legendStrings, nextLegendString];
        
        timestepIdx = timestepIdx + 1;
    end
    
    legend(legendStrings, 'Location', 'southeast')
    title(method);
    methodIdx=methodIdx+1;
end;

hold off;

% clf;
% plot(T, f(T));