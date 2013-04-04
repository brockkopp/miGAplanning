%% Initialize
% clc; clear all; 

%% World
inc = 0.01;     % cn resolution
xDim = 2/inc; 	%2m
yDim = 2/inc;   %2m
world = zeros(xDim,yDim);

armLen = [ 100 50 ];
armAng = [10 30;
          80 80];

armPos = updateArmPos(armLen, armAng(1,:));

obstacles = [ 50  80 100 120;
             120 150  30  50;
             140 160 140 180];

numObs = length(obstacles(:,1));
    
for obs=1:numObs
    for i=1:xDim
        for j=1:yDim
            if(i > obstacles(obs,1) && i < obstacles(obs,2) && j > obstacles(obs,3) && j < obstacles(obs,4))
                world(i,j) = 1;
            end
        end
    end
end



%% Calcular Configuration Space
% cSpaceLimits = [1 90 1 360];
% [cSpace F] = buildCspace(armLen, world, cSpaceLimits);


%% Plots
% plotWorld(world, armPos);
% plotCspace( cSpace );

plotAll( world, armPos, cSpace );

%% GA
[CxDim,~] = size(cSpace);
X = 1:CxDim;
% dna = [2 -0.1 0.1 0];
% [Y, fit] = splineFitness(dna,cSpace,armAng,X);

% FitnessFunction = @(xRow) splineFitness(xRow,cSpace,armAng,X);
% 
% % [Y, fit] = FitnessFunction(dna);
% % plot(X,Y);
% 
% gaoptions = gaoptimset('PopulationSize',100);
% 
% [dna, fval] = ga(FitnessFunction,4,[],[],[],[],[],[],[],gaoptions)

% offset = armAng(1,2) - (dna(4)*armAng(1,1)^4 + dna(3)*armAng(1,1)^3 + dna(2)*armAng(1,1)^2 + dna(1)*armAng(1,1));
% adj = offset*ones(size(X));
% Y = dna(4)*X.^4 + dna(3)*X.^3 + dna(2)*X.^2 + dna(1)*X + adj; 
% 
% plot(X,Y);
