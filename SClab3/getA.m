function [ A ] = getA( Nx, Ny )
%Construct A matrix

hx = 1/(Nx-1);
hy = 1/(Ny-1);

A = zeros(Nx*Ny);
for row = 1:size(A,1)
    i = mod((row-1),Nx)+1;
    j = fix((row-1)/Nx)+1;
    col = @(i,j) i+(j-1)*Nx;
    A(row, col(i,j)) = -2*(hx^2+hy^2);      %i,j
    if i>1 && j>1 && i<Nx && j<Ny           %Boundry condition
        A(row, col(i-1,j)) = hy^2;          %i-1,j
        A(row, col(i+1,j)) = hy^2;          %i+1,j
        A(row, col(i,j-1)) = hx^2;          %i-1,j
        A(row, col(i,j+1)) = hx^2;          %i+1,j
    end
end

end

