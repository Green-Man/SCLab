clear all;clf;hold on;
disp('==================')

%Time boundaries
t_start = 0;
dts = fliplr([1/1 1/2 1/4 1/8]);
t_end = 5;

%Given ODE and initial condition
dpdt = @(p) (1 - p/10)*p;
p0 = 1;

methods = {'Euler', 'Heun', 'Runge-Kutta'};
i = 1;
exactError = zeros(size(methods,2), size(dts,2));
approxError = zeros(size(methods,2), size(dts,2));

for method = methods
    %Plot each method on the separate subplot
    figure(1)
    hold on
    subplot(size(methods,2),1,i);
    legendStrings = {'Analytic'};
    
    %Analytical solution and its plot
    p = @(t) 10 ./ (1+9*exp(-t));
    P = p(t_start:dts(1):t_end);
    plot(t_start:dts(1):t_end, P, 'Color', [0.3 0.3 0.3]);
    
    nsBest = [];	%the most precise solution
    j=1;
    for dt = (dts)

        T = t_start:dt:t_end;
        
        %Execution of numerical methods
        ns = numerical( dpdt, p0, t_start, dt, t_end , method);
        
        %Save the most precise solution
        if dt == dts(1)
            nsBest = ns;
        end
        
        %Plots of numerical methods
        gColor = [0,0,0]; gColor(i) = min([dt^0.5 1]);
        hold on;
        plot(T, ns, '^', 'Color', gColor);
        hold on;
        
        interpNs = interp1(t_start:dt:t_end, ns, t_start:dts(1):t_end, 'linear');
        
        %Compute exact error
        exactErrorF = @(N) sqrt(sum((N-P).^2)*dt/5);
        exactError(i,j) = exactErrorF(interpNs);
        
        %Compute approximation error
        approximationErrorF = @(N) sqrt(sum((N-nsBest).^2)*dt/5);
        if dt ~= dts(1)
            approxError(i,j) = approximationErrorF(interpNs);
        end
        
        %Format a legend
        errorString = sprintf('exactError=%2.6f', exactError(i,j));
        dtString = sprintf('dt=%2.3f', dt);
        nextLegendString = sprintf('%s', dtString);
        legendStrings = [legendStrings, nextLegendString];
        
        j = j + 1;
    end
    
    legend(legendStrings, 'Location', 'southeast')
    title(method);
    i=i+1;
end

%Display exact error
disp('Exact error')
disp([exactError])


%compute reduction error when "decreasing dt at two times"
re = ones(size(exactError,1),size(exactError,1)-2);
for method = 1:size(exactError,1)
    for dt = 1:size(dts,2)-1
        re(method, dt) = exactError(method, dt+1)/exactError(method, dt);
    end
end
disp('Reduction error')
disp(re)
for i = 1:size(methods,2)
    subplot(size(methods,2),1,i);
    title([methods(i) sprintf('Mean error reduction %2.2f', mean(re(i,:),2))])
end

%Compute approximation error
disp('Approximation error')
disp(approxError)

% exactError(i,j) = exactErrorF(ns);





hold off;