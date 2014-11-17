 clc;
xmin = 0; ymin = 0;
xmax = 1; ymax = 1;
NJacobi = 32;
Nx = 20; hx = (xmax-xmin)/(Nx-1);
Ny = 20; hy = (ymax-ymin)/(Ny-1);

[X,Y] = meshgrid(xmin:hx:xmax,...
                 ymin:hy:ymax);
tAnalytic = @(x,y) sin(pi.*x).*sin(pi.*y);
% clf;figure(1);
% surf(X,Y,tAnalytic(X,Y));

%Construct B column vector
%TODO: construct without iteration
fB = @(row) sin(pi*mod((row-1),Nx)*hx)*sin(pi*fix((row-1)/Nx)*hy);
B = ones(Nx*Ny,1);
for row = 1:length(B)
    B(row) = fB(row);
end;
B = B*-2*pi^2*hx^2*hy^2;

Xt = jacobiSolver(B, Nx, Ny);

% Xt = A\B;
toc
Xt = reshape(Xt, [Nx Ny]);
Xt = Xt';

figure(2);clf;
mesh(X,Y,Xt);