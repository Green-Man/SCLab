 clc;
xmin = 0; ymin = 0;
xmax = 1; ymax = 1;
NJacobi = 32;
Nx = 20; hx = (xmax-xmin)/(Nx-1);
Ny = 20; hy = (ymax-ymin)/(Ny-1);
accuracy = 1e-4;

[X,Y] = meshgrid(xmin:hx:xmax,...
                 ymin:hy:ymax);
tAnalytic = @(x,y) sin(pi.*x).*sin(pi.*y);
% clf;figure(1);
% surf(X,Y,tAnalytic(X,Y));

%Construct A matrix
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

%Construct B column vector
%TODO: construct without iteration
fB = @(row) sin(pi*mod((row-1),Nx)*hx)*sin(pi*fix((row-1)/Nx)*hy);
B = ones(Nx*Ny,1);
for row = 1:length(B)
    B(row) = fB(row);
end;
B = B*-2*pi^2*hx^2*hy^2;

Xt = zeros(size(A,1),1);
tic
k = 1;
Anij = A - eye(size(A)).*A;
while 1
    previousNorm = norm(Xt);
    for i = 1:size(A, 1)
        Xt(i) = ( B(i) - sum(Anij(i,:)*Xt) )/A(i,i);
    end
    if abs(norm(Xt) - previousNorm) < accuracy
        break
    end;
    k = k + 1;
end
% Xt = A\B;
toc
Xt = reshape(Xt, [Nx Ny]);
Xt = Xt';

% figure(2);clf;
% mesh(X,Y,Xt);