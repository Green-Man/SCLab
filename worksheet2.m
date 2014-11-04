clc;clf;hold on;
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

t_end = 5*10^0*1;
t_end_plot = t_end;
newtonE = 1e-4;

%Given ODE and initial condition
dpdt = @(p) 7.*(1 - p./10).*p;
dpdt_prime = @(p) 7.*(1 - p./5);
p0 = 20;

%Analytical solution
p = @(t) 200 ./ (20-10*exp(-7*t));
exactT = t_start:dts(1):t_end;
P = p(exactT);

methods = {'Implicit Euler', 'Adams-Moulton'};
methods = {'Implicit Euler',...
           'Adams-Moulton',...
           'Adams-Moulton L1',...
           'Adams-Moulton L2'};
methodIdx = 1;
exactError = zeros(size(methods,2), size(dts,2));
approxError = zeros(size(methods,2), size(dts,2));
times = zeros(size(methods,2), size(dts,2));

for method = methods
    graphs = zeros(size(dts,2)+1, 1);
    
    %Plot each method on the separate subplot
    subplot(2,ceil(size(methods,2)/2),methodIdx);
    legendStrings = cell(1,length(dts)+1);
    legendStrings(1) = {'Analytic'};
    g = plot(exactT(1:round(t_end_plot/dts(1))),...
                    P(1:round(t_end_plot/dts(1))),...
                    'Color', [0.3 0.3 0.3]);
    graphs(1) = g(1);
    nsBest = [];	%the most precise solution
    timestepIdx=1;
    
    for dt = dts
        T = t_start:dt:t_end;
        hold on; ylim([0 20]);
        gColor = hsv2rgb([timestepIdx/length(dts) 0.95 0.6]);
        
        %Execution of numerical methods
        [ns,...
         times(methodIdx,timestepIdx),...
         e] = numerical(dpdt, dpdt_prime, p0,...
                        t_start, dt, t_end , method, newtonE);
        
        e = 1/realmin('double')*(1.-e)+e.*ns;
        
        %Save the most precise solution
        if dt == dts(1) nsBest = ns; end
        
        %Plot of the numerical methods
        plot(T(1:round(t_end_plot/dt)), e(1:round(t_end_plot/dt)),...
            '^',...
            'Color', gColor,...
            'MarkerSize',10);
        g = plot(T(1:round(t_end_plot/dt)), ns(1:round(t_end_plot/dt)),...
            '-gs',...
            'Color', gColor,...
            'LineWidth', 1.2,...
            'MarkerSize',2);
        graphs(timestepIdx+1) = g(1);
        
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
        dtString = sprintf('dt=1/%.0f', 1/dt);
        nextLegendString = sprintf('%s', dtString);
        legendStrings(timestepIdx + 1) = {nextLegendString};
        
        timestepIdx = timestepIdx + 1;
    end
    
    legend(graphs, legendStrings, 'Location', 'southeast')
    title(method);
    methodIdx=methodIdx+1;
end;

hold off;
