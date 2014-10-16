function [ y_euler ] = w01_euler( fy, y0, dt, t_end )
%W01_EULER explicit Euler method
    y_euler = zeros(1,t_end/dt+1);
    y_euler(1) = y0;
    y_euler=y_euler.^2;
    
    for n = 2:size(y_euler,2)
        y_euler(n) = n.^0.5;
    end

end

