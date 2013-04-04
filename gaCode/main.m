close all
%% Simulation Params
numSimulations = 5;

%% GA Configuration Params
PopulationSize = 100;

Generations = 100;
obstacleWeight = 5;
lengthWeightFactor = 0.2; 
lineResolution = 1; % The line is checked this often for collisions against the configurations space map

crossoverFraction = 0.5; % fraction of population that will reproduce
eliteCount = 1;

tournamentSize = PopulationSize / 2;
selectionFunction = {@selectionTournament, tournamentSize};

%mutationRate = 0.01; % Chance of a mutation for a particular vector entry of chosen individual to mutate
%TerminationConvergenceTolerance = 50;
NumGensAvg = 5; %Termination tolerance is averaged over this many generations

%% Robot Start&End Point
startPt = [10; 100];
endPt = [80; 300];

%% Setup Plot

% Load obstacle Grid
obsGrid = importdata('obsGrid.mat');

mapPlot = figure; 
hold on;
[xDim yDim] = size(obsGrid);
axis([0 xDim 0 yDim]);
image(100*(1-obsGrid)');
colormap(gray);
plot(startPt(1), startPt(2), 'g*', 'MarkerSize',6);
plot(endPt(1), endPt(2), 'g*', 'MarkerSize',6);

for i=0:numSimulations
%% Number of variables in chromosome
nvars = 5;

%% Coefficient (Gene) cosntraints
low = zeros(nvars,1);
upp = zeros(nvars,1);
range = zeros(2,nvars);

min = -20;
max = 20;
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

%% Call GA function
options = gaoptimset('PopInitRange',range);
options = gaoptimset(options,'PopulationSize',PopulationSize);
options = gaoptimset(options,'PopInitRange',PopulationInitializationRange);
options = gaoptimset(options,'Generations',Generations);
options = gaoptimset(options,'PlotFcns',{@gaplotbestf, @gaplotstopping});
%options = gaoptimset(options,'TolFun',TerminationConvergenceTolerance);
%options = gaoptimset(options,'StallGenLimit',NumGensAvg);
options = gaoptimset(options,'SelectionFcn',selectionFunction);
options = gaoptimset('MutationFcn', @mutationadaptfeasible);
options = gaoptimset('CrossoverFraction', crossoverFraction);
options = gaoptimset('EliteCount', eliteCount);

[x, Fval, exitFlag, Output] = ga(@(x) AKfitness(x,startPt, endPt, obstacleWeight, lengthWeightFactor, lineResolution),nvars,[],[],A_linEq,b_linEq,low,upp,[],[],options);

fprintf('Function Value = %g\n', Fval);
format short;
fprintf('y = %.3g + (%.3g)x + (%.3g)x^2 + (%.3g)x^3 + (%.3g)x^4 \n', x(1), x(2), x(3), x(4), x(5));


t = startPt(1):0.1:endPt(1);


figure(mapPlot);
plot(t,x(1)+x(2)*t + x(3)*t.^2 + x(4)*t.^3 + x(5)*t.^4, 'r');

%fprintf('Length = (%.3g)', minLength(x, startPt, endPt));
end