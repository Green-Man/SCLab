function [ y ] = w01_huen( fy, y0, t0, dt, t_end )
%W01_HUEN method Heun
    y = zeros(1,size(0:dt:t_end,2));
    y(1) = y0;
    
    for n = 2:size(y,2)
        dydt1 = fy(y(n-1));
        dydt2 = fy(y(n-1)+dt*dydt1);
        y(n) = y(n-1) + dt*(dydt1+dydt2)/2;
    end


end

