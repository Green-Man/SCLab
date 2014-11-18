clc;
xmin = 0; ymin = 0;
xmax = 1; ymax = 1;
NJacobi = 32;
N = [7 15 31 63];
N = [7 15];

for ni = 1:length(N)
    Nx = N(ni); hx = (xmax-xmin)/(Nx-1);
    Ny = N(ni); hy = (ymax-ymin)/(Ny-1);

    [X,Y] = meshgrid(xmin:hx:xmax,...
                     ymin:hy:ymax);
    tAnalytic = @(x,y) sin(pi.*x).*sin(pi.*y);
    figure(1);
    subplot(2,2,ni); zlim([0,1]);
    surf(X,Y,tAnalytic(X,Y));

    %Construct B column vector
    %TODO: construct without iteration
    fB = @(row) sin(pi*mod((row-1),Nx)*hx)*sin(pi*fix((row-1)/Nx)*hy);
    B = ones(Nx*Ny,1);
    for row = 1:length(B)
        B(row) = fB(row);
    end;
    B = B*-2*pi^2*hx^2*hy^2;
    tic
%     Xt = getA(Nx, Ny, 'full') \ B; %Full matrix and direct solver
%     Xt = getA(Nx, Ny, 'sparse') \ B; %Sparse matrix and direct solver
    Xt = gaussSeidelSparse(B, Nx, Ny); % Achtung! Takes Around one hour to compute 63x63 mesh
    toc
    
    Xt = reshape(Xt, [Nx Ny]);
    Xt = Xt';

    figure(2);
    colormap('jet');
    sp = subplot(2,2,ni);
    mesh(X,Y,Xt);
    zlim(sp,[0 1]);
    
    figure(3);
    subplot(2,2,ni);
    v = linspace(0,1,6);
    [C,h] = contour3(X,Y,Xt, v);
    clabel(C,h);

end