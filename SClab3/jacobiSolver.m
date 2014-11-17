function [ X ] = jacobiSolver( B, Nx, Ny )
%JACOBISOLVER Summary of this function goes here
%   Detailed explanation goes here

accuracy = 1e-4;
A = getA(Nx, Ny);

X = zeros(size(A,1),1);
tic
k = 1;
Anij = A - eye(size(A)).*A;
while 1
    previousNorm = norm(X);
    for i = 1:size(A, 1)
        X(i) = ( B(i) - sum(Anij(i,:)*X) )/A(i,i);
    end
    if abs(norm(X) - previousNorm) < accuracy
        break
    end;
    k = k + 1;
end

end

