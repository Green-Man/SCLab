function [ Z ] = implicitEuler( Nx, Ny, dt,  Zp )
    
    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    accuracy = 1e-6;
    R = 1;
    d = 1 + 2*dt*(1/hx^2+1/hy^2);
    Z = Zp;
    while abs(R)>accuracy
        for i = 2:Nx+1
            for j = 2:Ny+1
                d2 = (Z(i-1,j)+Z(i+1,j))/hx^2+...
                         (Z(i,j-1)+Z(i,j+1))/hy^2;
                Z(i,j) = (Zp(i,j)+dt*d2)/d;
            end
        end
        
        %Calculating residual
        R = 0;
        for i = 2:Nx+1
            for j = 2:Ny+1
                r = (Z(i-1,j)+Z(i+1,j))/hx^2+...
                         (Z(i,j-1)+Z(i,j+1))/hy^2-...
                         2*Z(i,j)*(1/hx^2+1/hy^2);
                r = r*dt + Zp(i,j) - Z(i,j);
                R = R + r^2;
            end
        end
        R = sqrt( R/(Nx*Ny));
    end
end

