clear all;clf;hold on;
disp('============================================')

%Time boundaries
t_start = 0;
dts = fliplr([1/1 1/2 1/4 1/8]); % Flip list to compute the most accurate solution first
t_end = 5;

%Given ODE and initial condition
dpdt = @(p) (1 - p/10)*p;
p0 = 1;

methods = {'Euler', 'Heun', 'Runge-Kutta'};
methodIdx = 1;
exactError = zeros(size(methods,2), size(dts,2));
approxError = zeros(size(methods,2), size(dts,2));

for method = methods
    %Plot each method on the separate subplot
    figure(1)
    hold on
    subplot(size(methods,2),1,methodIdx);
    legendStrings = {'Analytic'};
    
    %Analytical solution and its plot
    p = @(t) 10 ./ (1+9*exp(-t));
    exactT = t_start:dts(1):t_end;
    P = p(exactT);
    plot(exactT, P, 'Color', [0.3 0.3 0.3]);
    
    nsBest = [];	%the most precise solution
    timestepIdx=1;
    for dt = dts
        
        T = t_start:dt:t_end;
        
        %Execution of numerical methods
        ns = numerical( dpdt, p0, t_start, dt, t_end , method);
        
        %Save the most precise solution
        if dt == dts(1)
            nsBest = ns;
        end
        
        %Plots of numerical methods
        gColor = [0,0,0]; gColor(methodIdx) = min([dt^0.5 1]);
        hold on;
        plot(T, ns, '^', 'Color', gColor);
        hold on;

        %Compute exact error
        interpP = interp1(exactT, P, t_start:dt:t_end, 'linear');
        exactErrorF = @(N) sqrt(sum((N-interpP).^2)*dt/5);
        exactError(methodIdx,timestepIdx) = exactErrorF(ns);
        
        %Compute approximation error
        interpNsBest = interp1(exactT, nsBest, t_start:dt:t_end, 'linear');
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
end

%Display exact error
disp('Exact error')
disp(exactError)


%compute reduction error when "decreasing dt two times"
re = ones(size(exactError));
for method = 1:size(exactError,1)
    for dtIdx = 2:size(exactError,2)
        re(method, dtIdx) = exactError(method, dtIdx)/exactError(method, dtIdx-1);
    end
end
disp('Errors reduction')
format shortG
disp(re)

for methodIdx = 1:size(methods,2)
    subplot(size(methods,2),1,methodIdx);
    title([methods(methodIdx) sprintf('Mean error reduction %2.2f', mean(re(methodIdx,:),2))])
end

%Display table of approximation errors
disp('Approximation errors')
disp(approxError)

hold off;