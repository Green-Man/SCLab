function [ y_euler ] = w01_euler( fy, y0, t0, dt, t_end )
%W01_EULER explicit Euler method
    y_euler = zeros(1,size(0:dt:t_end,2));
    y_euler(1) = y0;
    
    for n = 2:size(y_euler,2)
        dpdt = fy(y_euler(n-1));
        y_euler(n) = y_euler(n-1) + dt*dpdt;
    end

end

