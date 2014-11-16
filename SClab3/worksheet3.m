clf; clc;
hold on;
xmin = 0; ymin = 0;
xmax = 1; ymax = 1;
NJacobi = 32;
Nx = 2; hx = (xmax-xmin)/Nx;
Ny = 3; hy = (ymax-ymin)/Ny;

[X,Y] = meshgrid(xmin:hx:xmax,...
                 ymin:hy:ymax);
tAnalytic = @(x,y) sin(pi.*x).*sin(pi.*y);
figure(1);
surf(X,Y,tAnalytic(X,Y));

%Construct A matrix
A = zeros((Nx+1)*(Ny+1));
for row = 1:size(A,1)
    row
    i = fix((row-1)/(Nx+1))+1;
    disp(i)
    j = mod((row-1),(Ny+1))+1;
    disp(j)
    disp((i-1)*(Nx+1)+(j-1)+1)
    A(row, (i-1)*(Nx+1)+(j-1)+1) = -2*(hx^2+hy^2);      %i,j
    if i>1 && j>1 && i<(Nx+1) && j<(Ny+1)               %Boundry condition
        A(row, (i-2)*(Nx+1)+(j-1)+1) = hy^2;            %i-1,j
        A(row, (i)*(Nx+1)+(j-1)+1) = hy^2;              %i+1,j
        A(row, (i-1)*(Nx+1)+(j-2)+1) = hx^2;            %i-1,j
        A(row, (i-1)*(Nx+1)+(j)+1) = hx^2;              %i+1,j
    end
end

%Construct B column vector
%TODO: construct without iteration
fB = @(x,y) sin(pi*x)*sin(pi*y);
B = ones((Nx+1)*(Ny+1),1);
for i = 1:length(B)
    B(i) = fB(fix((i-1)/(Nx+1))*hx, mod((i-1),(Ny+1))*hy);
end;
B = B*-2*pi^2*hx^2*hy^2;

Xt = zeros(size(A,1),1);
for k = 1:NJacobi
    for i = 1:size(A, 1)
        s = 0;
        for j=1:size(A,1)
            if i~=j
                s = s + A(i,j)*Xt(j);
            end;
        end;
        Xt(i) = ( B(i) - s' )/A(i,i);
    end
end
Xt = A\B;
Xt = reshape(Xt, [Nx+1  Ny+1]);
Xt = Xt';
figure(2);clf;
mesh(X,Y,Xt);