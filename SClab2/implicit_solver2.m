%Implicit solver 


function Y = implicit_solver2(p0,dt,T,eps_min,nmax,fun_newton,fun_newton_prim,p_fun,method)

    % p0 - initial value,
    % dt - timestep,
    % T - T end,
    % eps_min - border error in Newton's Method,
    % nmax - maximal number of steps for Newton's method,
    % fun_newton - Function needed for Newton's method,
    % fun_newton_prim - Derivative of the above function "p_fun_newton" also 
    %                   needed for Newton's method,
    % p_fun - basic equation (y' = f(y)),
    % method - name of used method ('euler' for Impliciit Euler,
    %                               'adams' for Adam Mouldon Method).

    t = 0:dt:T; %time
    newton_0 = 20; %Start point of Newton's method
    
switch (method),
    
    
    case char('euler'),  
        p_euler = zeros(1,length(t));
        p_euler(1) = p0;

        for i=1:length(t)-1,
            [p_solution, if_found] = newton_solver(newton_0,p_euler(i),dt,eps_min,nmax,fun_newton,fun_newton_prim);
            
            %if_found determines whether it existed solution for specific polynomial
            %if not (= 0), then whole solution is erased and we break out
            if if_found == 0, 
                p_euler = NaN(1,length(t));
                break;
            end   
            
            p_euler(i+1) = p_solution;
            newton_0 = p_solution;
        end
        Y = p_euler;
        
    case char('adams'),
        p_adams = zeros(1,length(t));
        p_adams(1) = p0;

        for i=1:length(t)-1,
            [p_solution, if_found] = newton_solver(newton_0,p_adams(i),dt,eps_min,nmax,fun_newton,fun_newton_prim);
            
            %if_found determines whether it existed solution for specific polynomial
            %if not (= 0), then whole solution is erased and we break out
            if if_found == 0,
                p_adams = NaN(1,length(t));
                break;
            end  
            
            p_adams(i+1) = p_solution;
            newton_0 = p_solution;
        end
        Y = p_adams;
      
end
        
end