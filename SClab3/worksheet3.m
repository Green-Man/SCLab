clf; clc;
hold on;
xmin = 0; ymin = 0;
xmax = 1; ymax = 1;
step = 0.01;
NJacobi = 30;
Nx = 3; dx = (xmax-xmin)/Nx;
Ny = 3; dy = (ymax-ymin)/Ny;

[X,Y] = meshgrid(xmin:dx:xmax,...
                 ymin:dy:ymax);
tAnalytic = @(x,y) sin(pi.*x).*sin(pi.*y);
% surf(X,Y,tAnalytic(X,Y));

%Construct A matrix

%Construct B column vector
fB = @(x,y) -2*pi^2*sin(pi*x)*sin(pi*y);
B = ones((Nx+1)*(Ny+1),1);
for i = 1:length(B)
    B(i) = fB(fix((i-1)/(Nx+1))*dx, mod((i-1),(Ny+1))*dy);
end;

X = zeros(size(A,1),1);
for k = 1:NJacobi
    for i = 1:size(A, 1)
        s = 0;
        for j=1:size(A,2)
            if i~=j
                s = s + A(i,j)*X(j);
            end;
        end;
        X(i) = ( B(i) - s' )/A(i,i);
    end
end