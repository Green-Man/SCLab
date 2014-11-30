clc;

plotTimes = [0 1/8 2/8 3/8 4/8]; endTime = plotTimes(end);
dTs = 1./[64 128 256 512 1024 2048 4096];
spatialRes = [3 7 15 31];
dims = [spatialRes;spatialRes];

% tables init
rownames = {'runtime', 'storage'}; rownames2 = {'error', 'error red.'};
varnames = {}; varnames2 = {};
data1 = []; data2 = []; data3 = []; data4 = [];
meshIdx = 0;
for dim=dims
    Nx = dim(1,:); Ny = dim(2,:);
    hx = 1/(Nx+1);hy = 1/(Ny+1);
    [XX,YY] = meshgrid(0:hx:1, 0:hy:1);
    dti = 1;
    %Explicit Euler
    for dt = dTs
        %Initial condition
        t = 0;
        Zp = ones(Nx+2, Ny+2);
        Zp(:,1) = 0;Zp(:,end) = 0;
        Zp(1,:) = 0;Zp(end,:) = 0;
        while t <= endTime
            if t == 0
                Z = Zp;
            else
                Z = explicitEuler(Nx, Ny, dt, Zp);
            end
            Zp = Z;

            %plotting EE
            [isPlottingTime, fi] = ismember(t, plotTimes);
            if isPlottingTime
                f = figure(fi);
                f.Name = sprintf('t = %1.0f/8', t*8);
                subplot(length(spatialRes),...
                        length(dTs),...
                        meshIdx*length(dTs) + dti);
                mesh(XX,YY, Z);
                title(sprintf('dt=1/%4.0f, Nxy=%i',1/dt,spatialRes(meshIdx+1)));
                axis([0 1 0 1 0 max([1e-3 max(Z(:))])]);
            end
            t = t + dt;
        end
        dti = dti + 1;
    end
    
    %Implicit Euler
    dt =1/64;
    t = 0;
    Zp = ones(Nx+2, Ny+2);
    Zp(:,1) = 0;Zp(:,end) = 0;
    Zp(1,:) = 0;Zp(end,:) = 0;

    while t <= endTime
        if t == 0
            Z = Zp;
        else
            Z = implicitEuler(Nx, Ny, dt, Zp);
        end
        Zp = Z;
        
        %Plot
        [isPlottingTime, fi] = ismember(t, plotTimes);
        if isPlottingTime
            f = figure(length(plotTimes)+1);
            f.Name = sprintf('Implicit Euler, dt = 1/64');
            subplot(length(spatialRes),...
                    length(plotTimes),...
                    meshIdx*length(plotTimes) + fi);
            mesh(XX,YY, Z);
            title(sprintf('dt = 1/%4.0f', 1/dt));
            axis([0 1 0 1 0 max([1e-3 max(Z(:))])]);
        end
        t = t + dt;
    end
    meshIdx = meshIdx + 1;
end

% tables output
% disp('iterative solution with Gauss-Seidel'); disp(array2table(data3, 'RowNames', rownames, 'VariableNames', varnames2));
