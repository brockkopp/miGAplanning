% reference http://www.mathworks.com/help/gads/genetic-algorithm-options.html

close all
%% Simulation Params
numSimulations = 3;

%% GA Configuration Params
PopulationSize = 125;
Generations = 25;

coefRangeMin = -100;
coefRangeMax = 100;

obstacleWeight = 2;
lengthWeightFactor = 0.01; 
lineResolution = 1; % The line is checked this often for collisions against the configurations space map
jerkWeight = 0*0.0001;

TerminationConvergenceTolerance = 0.1;
NumGensAvg = 10;

crossoverFraction = 0.80; % fraction of population that will reproduce
eliteCount = 1; %floor(PopulationSize * 0.01);
%  ratio = 0.5;
% crossoverFunction = {@crossoverintermediate, ratio};
crossoverFunction = @crossoverheuristic;
% crossoverFunction = @crossovertwopoint;
mutationFunction = @mutationadaptfeasible;

fitnessScalingFunction = @fitscalingprop;

% tournamentSize = PopulationSize * 0.1;
tournamentSize = 6;
selectionFunction = {@selectiontournament, tournamentSize};
%selectionFunction = @selectionstochunif;

%% Robot Start&End Point
startPt = [20; 50];
endPt = [80; 300];

%% Setup Plot

% Load obstacle Grid
% obsGrid = importdata('obsGrid.mat');
obsGrid = obsGrid; 

mapPlot = figure; 
hold on;
[xDim yDim] = size(obsGrid);
axis([0 xDim 0 yDim]);
image(100*(1-obsGrid)');
colormap(gray);
plot(startPt(1), startPt(2), 'g*', 'MarkerSize',6);
plot(endPt(1), endPt(2), 'g*', 'MarkerSize',6);

t = startPt(1):0.1:endPt(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Iterate through all simulations%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for j=1:numSimulations
%% Number of variables in chromosome
nvars = 5;

%% Coefficient (Gene) cosntraints
low = zeros(nvars,1);
upp = zeros(nvars,1);
range = zeros(2,nvars);

min = coefRangeMin;
max = coefRangeMax;

for(i = 1:nvars)
    low(i) = min;
    upp(i) = max;
    range(1,i) = min;
    range(2,i) = max;
end
PopulationInitializationRange = range;
%% Linear Equalities
x1 = startPt(1);
x2 = endPt(1);

y1 = startPt(2);
y2 = endPt(2);

A_linEq = zeros(2, nvars);
for(i = 1:nvars)
    A_linEq(1,i) = x1^(i-1);
    A_linEq(2,i) = x2^(i-1);
end
    
b_linEq = [y1; y2];

%% Generate initial population
% initialPopulation = coefRangeMin + (coefRangeMax - coefRangeMin)*rand(PopulationSize,nvars);
% initialPopulation(1,:) = [0 ((startPt(2) - endPt(2)) / (startPt(1) - endPt(1))) 0 0 0];  

%% Call GA function
options = gaoptimset('PopInitRange',range);
options = gaoptimset(options,'PopulationSize',PopulationSize);
options = gaoptimset(options,'PopInitRange',PopulationInitializationRange);
% options = gaoptimset(options,'InitialPopulation', initialPopulation);

options = gaoptimset(options,'Generations',Generations);
% options = gaoptimset(options,'PlotFcns',{@gaplotbestf, @gaplotstopping});
options = gaoptimset(options,'TolFun',TerminationConvergenceTolerance);
options = gaoptimset(options,'StallGenLimit',NumGensAvg);
options = gaoptimset(options,'SelectionFcn',selectionFunction);
options = gaoptimset(options,'MutationFcn', mutationFunction);
options = gaoptimset(options,'CrossoverFraction', crossoverFraction);
options = gaoptimset(options,'CrossoverFcn', crossoverFunction);
options = gaoptimset(options,'EliteCount', eliteCount);
options = gaoptimset(options,'FitnessScalingFcn', fitnessScalingFunction);

[x, Fval, exitFlag, Output] = ga(@(x) AKfitness(x,startPt, endPt, obstacleWeight, lengthWeightFactor, jerkWeight, lineResolution, j),nvars,[],[],A_linEq,b_linEq,low,upp,[],[],options);

fprintf('Fitness Value = %g\n', Fval);
fprintf('Generations   = %g\n', Output.generations);
format short;
fprintf('y = %.3g + (%.3g)x + (%.3g)x^2 + (%.3g)x^3 + (%.3g)x^4 \n', x(1), x(2), x(3), x(4), x(5));

A = x(1);
B = x(2);
C = x(3);
D = x(4);
E = x(5);

collision = 1;
offScreen = 0;
numCollisions = 0;
djerk = 0;
ddist = 0;
dcoll = 0;

y2 = @(t) sqrt(1 + (B + 2*C*t + 3*D*t.^2 + 4*E*t.^3).^2); % (B^2 + 1) + (4*B*C)*t + (4*C^2)*t.^2);
ddist = integral(y2, startPt(1), endPt(1)) * lengthWeightFactor;
for i=startPt(1):lineResolution:endPt(1)
    y = A + B*i + C*i^2 + D*i^3 + E*i^4;
    if (y > yDim || y < 0)
%         offScreen = 1;
        numCollisions = 10000;
        break;
    elseif (obsGrid(i,ceil(y)) == 1) %if within obstacle
        numCollisions = numCollisions + collision;
        dcoll = dcoll + obstacleWeight * ceil(abs(y_v));
    end
    y_v =   B   + 2*C*i + 3*D*i^2 + 4*E*i^3;
%     y_a =         2*C   + 6*D*i   + 12*E*i^2;
    y_jerk =                 6*D     + 24*E*i;
    djerk = djerk + y_jerk*jerkWeight;
    
end
disp 'jerk'
djerk

disp 'Obstacle'
dcoll

disp 'length'
ddist

fprintf('Number of Collisions: %g\n', numCollisions);

figure(mapPlot);
plot(t,x(1)+x(2)*t + x(3)*t.^2 + x(4)*t.^3 + x(5)*t.^4, 'r');

%fprintf('Length = (%.3g)', minLength(x, startPt, endPt));

end

disp '_Done'