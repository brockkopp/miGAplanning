%% Initialize
% clc; clear all; 

%% World
inc = 0.01;     % cn resolution
xDim = 2/inc; 	%2m
yDim = 2/inc;   %2m
world = zeros(xDim,yDim);

numObs = 5;
obsSize = [10 20];

armBase = [ 100 100];
armLen = [50 50];
armAng = [10 30];

armPos = updateArmPos(armBase, armLen, armAng(:));
obstacles = zeros(numObs,4);

for i=1:numObs
    rad = round(obsSize(1) + (obsSize(2)-obsSize(1))*rand(1,1));
    origX = round(10 + (xDim-20)*rand(1,1));
    origY = round(10 + (yDim-20)*rand(1,1));
    
    obstacles(i,:) = [ origX-rad origX+rad origY-rad origY+rad ];
end

for obs=1:numObs
    for i=1:xDim
        for j=1:yDim
            if(i > obstacles(obs,1) && i < obstacles(obs,2) && j > obstacles(obs,3) && j < obstacles(obs,4))
                world(i,j) = 1;
            end
        end
    end
end

% plotWorld(world, armPos);

%% Calcular Configuration Space
% if(~exist('cSpace','var'))
    cSpaceLimits = [1 720 1 360];
    [cSpace F] = buildCspace(armBase, armLen, world, cSpaceLimits);
% end

%% Plots
% plotWorld(world, armPos);
% plotCspace( cSpace );
plotAll( world, armPos, cSpace );
