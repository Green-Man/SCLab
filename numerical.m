function [ y, t, e] = numerical( fy, fy_prime, y0,...
                                t_start, dt, t_end,...
                                method, newtonE )
%NUMERICAL methods
    tic
    T = [t_start:dt:t_end];
    
    y = zeros(1,size(T,2));
    e = y; % Row to indicate if an error occured
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
            y(n) = nsolve(f, fp, y(n-1), newtonE);
        end
        
    elseif  strcmp(method,'Adams-Moulton')
        for n = 2:size(y,2)
            f =  @(x) y(n-1)+dt/2*( fy(x) + fy(y(n-1)) ) - x;
            fp = @(x) dt/2*fy_prime(x) - 1;
            %Define the root(s) existance condition
            D = -49*dt^2*y(n-1)*y(n-1) + 490*dt^2*y(n-1)+...
                    1225*dt^2+140*dt*y(n-1)-700*dt+100 >= 0;
            if D 
                %Regular Newton's method solution
                y(n) = nsolve(f, fp, y(n-1), newtonE);
            else
                %Adams-Moulton L1
                y(n) = y(n-1)*(20-7*dt*(y(n-1)-20));
                y(n) = y(n)/(7*dt*y(n-1)+20);
                e(n) = 1;
            end
        end
            
    elseif  strcmp(method,'Adams-Moulton L1')
        for n = 2:size(y,2)
            y(n) = y(n-1)*(20-7*dt*(y(n-1)-20));
            y(n) = y(n)/(7*dt*y(n-1)+20);
        end
    
    elseif  strcmp(method,'Adams-Moulton L2')
        for n = 2:size(y,2)
            y(n) = y(n-1)*(20-7*dt*(y(n-1)-10));
            y(n) = y(n)/(7*dt*(y(n-1)-10)+20);
        end
    end
    t = toc;
end

