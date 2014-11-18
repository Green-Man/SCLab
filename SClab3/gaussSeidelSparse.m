function [ X ] = gaussSeidelSparse( B, Nx, Ny )
% Solve without storing the system matrix
    
    accuracy = 1e-4;
    X = zeros(1, Nx*Ny);
    i = 1;
    while 1
        previousNorm = norm(X);
        for row = 1:Nx*Ny
            Arow = getRow(row, Nx, Ny);
            s1 = dot(Arow(1:row-1), X(1:row-1));
            s2 = dot(Arow(row+1:end), X(row+1:end));
            X(row) = ( B(row) - s1 -s2 )/Arow(row);
        end
        if abs(norm(X) - previousNorm) < accuracy
            break
        end;
        i = i + 1;
    end
    disp(i);
end

