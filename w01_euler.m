function [ y ] = w01_euler( fy, y0, t0, dt, t_end )
%W01_EULER explicit Euler method
    y = zeros(1,size(0:dt:t_end,2));
    y(1) = y0;
    
    for n = 2:size(y,2)
        dydt = fy(y(n-1));
        y(n) = y(n-1) + dt*dydt;
    end

end

