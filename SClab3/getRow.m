function [ A ] = getRow( row, Nx, Ny )
%Construct one row of A matrix

    hx = 1/(Nx-1);
    hy = 1/(Ny-1);

    A = zeros(1, Nx*Ny);

    i = mod((row-1),Nx)+1;
    j = fix((row-1)/Nx)+1;
    col = @(i,j) i+(j-1)*Nx;
    A(col(i,j)) = -2*(hx^2+hy^2);      %i,j
    if i>1 && j>1 && i<Nx && j<Ny      %Boundry condition
        A(col(i-1,j)) = hy^2;          %i-1,j
        A(col(i+1,j)) = hy^2;          %i+1,j
        A(col(i,j-1)) = hx^2;          %i-1,j
        A(col(i,j+1)) = hx^2;          %i+1,j
    end


end

