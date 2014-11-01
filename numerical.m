function [ y t ] = numerical( fy, fy_prime, y0, t_start, dt, t_end, method, newtonE )
%NUMERICAL methods
    tic
    y = zeros(1,size(t_start:dt:t_end,2));
    y(1) = y0;
    
    if strcmp(method,'Runge-Kutta')
        for n = 2:size(y,2)
            dydt1 = fy( y(n-1) );
            dydt2 = fy( y(n-1)+dt/2*dydt1 );
            dydt3 = fy( y(n-1)+dt/2*dydt2 );
            dydt4 = fy( y(n-1)+dt*dydt3 );
            y(n) = y(n-1) + dt*(dydt1+2*dydt2+2*dydt3+dydt4)/6;
        end
        
	elseif  strcmp(method, 'Heun')
        for n = 2:size(y,2)
            dydt1 = fy(y(n-1));
            dydt2 = fy(y(n-1)+dt*dydt1);
            y(n) = y(n-1) + dt*(dydt1+dydt2)/2;
        end
        
	elseif  strcmp(method,'Explicit Euler')
        for n = 2:size(y,2)
            dydt = fy(y(n-1));
            y(n) = y(n-1) + dt*dydt;
        end
    
 	elseif  strcmp(method,'Implicit Euler')
        for n = 2:size(y,2)
            
            f =  @(x) y(n-1)+fy(x)*dt-x;
            fp = @(x) fy_prime(x)*dt-1;
            y_predict = nsolve(f, fp, newtonE);
            y(n) = y(n-1) + dt*fy(y_predict);
        end
    end
    t = toc;
end

