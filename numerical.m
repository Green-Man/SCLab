function [ y ] = numerical( fy, y0, t0, dt, t_end, method )
%NUMERICAL methods
    y = zeros(1,size(0:dt:t_end,2));
    y(1) = y0;
    
    if strcmp(method,'rk')
        for n = 2:size(y,2)
            dydt1 = fy( y(n-1) );
            dydt2 = fy( y(n-1)+dt/2*dydt1 );
            dydt3 = fy( y(n-1)+dt/2*dydt2 );
            dydt4 = fy( y(n-1)+dt*dydt3 );
            y(n) = y(n-1) + dt*(dydt1+2*dydt2+2*dydt3+dydt4)/6;
        end
        
	elseif  strcmp(method, 'heun')
        for n = 2:size(y,2)
            dydt1 = fy(y(n-1));
            dydt2 = fy(y(n-1)+dt*dydt1);
            y(n) = y(n-1) + dt*(dydt1+dydt2)/2;
        end
        
	elseif  strcmp(method,'euler')
        for n = 2:size(y,2)
            dydt = fy(y(n-1));
            y(n) = y(n-1) + dt*dydt;
        end
end

