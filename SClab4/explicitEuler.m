function [ Z ] = explicitEuler( Nx, Ny, dt,  Zp )

    hx = 1/(Nx+1);
    hy = 1/(Ny+1);
    D = zeros(Nx+2,Ny+2);
    for i = 2:Nx+1
        for j = 2:Ny+1
            D(i,j) = (Zp(i-1,j)-2*Zp(i,j)+Zp(i+1,j))/hx^2+...
                     (Zp(i,j-1)-2*Zp(i,j)+Zp(i,j+1))/hy^2;
        end
    end
    Z = Zp+dt*D;

end

