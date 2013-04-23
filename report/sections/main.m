%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     MATLAB CODE - main.m     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Purpose: Executes genetic algorim to determine path across pre-defined envrionment. Used to collect datasets.

% Initialize file variables
simName = 'longtest';
simdir = strcat('sim_', simName);
mkdir(simdir);
fileName = strcat('gaData', int2str(runNumber),'_',datestr(now,'dd.mm.yyyy-HH.MM.SS'));

%% Initialize simulation params
numCspaces = 3;
numPointSets = 5;
colors = ['b' 'r' 'g' 'k' 'c'];
numSimulations = 500; %per points per cSpace
outputData = zeros(numSimulations, 17);

%% GA Configuration Params
PopulationSize = 75;
Generations = 50;

coefRangeMin = -500;
coefRangeMax = 500;

% GA Penalties (penalty weight relative)
obstacleWeight = 2;
lengthWeightFactor = 0.01;
lineResolution = 1;
jerkWeight = 1*0;

% GA Termination conditions
TerminationConvergenceTolerance = 0.001;
NumGensAvg = 10;

% GA Crossover parameters
crossoverFraction = 0.80;
eliteCount = 1;
crossoverFunction = @crossoverheuristic;
mutationFunction = @mutationadaptfeasible;
fitnessScalingFunction = @fitscalingprop;

% GA Crossover selection criteria
tournamentSize = 2;
selectionFunction = {@selectiontournament, tournamentSize};

% Setup Plot
cSpaceFilenames = [ 'cSpace2',
                    'cSpace3',
                    'cSpace4',
                    'obsGrid'];

pointSetVector = [  7, 104; 349, 112;
                   15, 142; 238, 247;
                    7, 263; 311, 159;
                   91, 222; 314,  24;
                   58, 219; 184, 341;
                  155, 318; 352, 146;
                  116, 194; 206, 131;
                  156, 259; 286,  83;
                   48, 257; 208, 112;
                  151, 248; 308, 296;
                  105,  91; 319, 107;
                  143, 328; 264, 41;
                  120, 135; 263, 3;
                  65 , 208; 359, 314;
                  115, 229; 336, 47;]

% Iterate through each environment
for cSpaceIteration = 1:numCspaces

    cSpaceID = cSpaceIteration;
    fpath = strcat('../configSpace/',  cSpaceFilenames(cSpaceIteration,:), '.mat')
    obsGrid = importdata(fpath);

% Plot environment
    mapPlot = figure;
    hold on;
    [xDim yDim] = size(obsGrid);
    axis([0 xDim 0 yDim]);
    image(100*(1-obsGrid)');
    colormap(gray);

% Iterate through each start/end point configuration
    for pointsIteration = 1:numPointSets
    pointSetID = pointsIteration;
    cSpace = obsGrid;
    startPt = [pointSetVector(pointsIteration*2-1 + (cSpaceIteration-1)*10, 1), pointSetVector(pointsIteration*2-1 + (cSpaceIteration-1)*10, 2)];
    endPt =   [pointSetVector(pointsIteration*2 + (cSpaceIteration-1)*10, 1), pointSetVector(pointsIteration*2 + (cSpaceIteration-1)*10, 2)];

% Plot start/end points
    figure(mapPlot);
    plot(startPt(1), startPt(2), 'g*', 'MarkerSize',6);
    plot(endPt(1), endPt(2), 'g*', 'MarkerSize',6);

    t = startPt(1):0.1:endPt(1);

% Repeat algorithm execution for dataset
        for j=1:numSimulations

% Number of variables in chromosome
        nvars = 4;

% Coefficient (Gene) cosntraints
        low = zeros(nvars,1);
        upp = zeros(nvars,1);
        range = zeros(2,nvars);

        mmin = coefRangeMin;
        mmax = coefRangeMax;

        for(i = 1:nvars)
            low(i) = mmin;
            upp(i) = mmax;
            range(1,i) = mmin;
            range(2,i) = mmax;
        end
        PopulationInitializationRange = range;

% Linear Equalities
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

% Define GA options
        options = gaoptimset('PopInitRange',range);
        options = gaoptimset(options,'PopulationSize',PopulationSize);
        options = gaoptimset(options,'PopInitRange',PopulationInitializationRange);
        options = gaoptimset(options,'Generations',Generations);
        options = gaoptimset(options,'TolFun',TerminationConvergenceTolerance);
        options = gaoptimset(options,'StallGenLimit',NumGensAvg);
        options = gaoptimset(options,'SelectionFcn',selectionFunction);
        options = gaoptimset(options,'MutationFcn', mutationFunction);
        options = gaoptimset(options,'CrossoverFraction', crossoverFraction);
        options = gaoptimset(options,'CrossoverFcn', crossoverFunction);
        options = gaoptimset(options,'EliteCount', eliteCount);
        options = gaoptimset(options,'FitnessScalingFcn', fitnessScalingFunction);
		
% Execute GA
        tic;
        [x, Fval, exitFlag, Output] = ga(@(x) AKfitness(x,startPt, endPt, obstacleWeight, lengthWeightFactor, jerkWeight, lineResolution, j),nvars,[],[],A_linEq,b_linEq,low,upp,[],[],options);
        gaLengthTime = toc;
		
% Print GA results
        fprintf('Fitness Value = %g\n', Fval);
        fprintf('Generations   = %g\n', Output.generations);
        format short;

        A = x(1);
        B = x(2);
        C = x(3);
        D = x(4);

        collision = 1;
        offScreen = 0;
        numCollisions = 0;
        djerk = 0;
        ddist = 0;
        dcoll = 0;
		
% Evaluate path
        y2 = @(t) sqrt(1 + (B + 2*C*t + 3*D*t.^2).^2);
        ddist = integral(y2, startPt(1), endPt(1)) * lengthWeightFactor;
        solutionLength = ddist;
        for i=startPt(1):lineResolution:endPt(1)
            y = A + B*i + C*i^2 + D*i^3;
            if (y > yDim || y < 0)
                numCollisions = 10000;
                dcoll = 10000;
                break;
            elseif (obsGrid(i,ceil(y)) == 1) %if within obstacle
                numCollisions = numCollisions + collision;
                dcoll = dcoll + obstacleWeight * ceil(abs(y_v));
            end
            y_v =   B   + 2*C*i + 3*D*i^2;
            y_jerk =                 6*D ;
            if (y_jerk > djerk)
                djerk = y_jerk;
            end
        end

        maxJerk = djerk;

% Plot trajectory
        figure(mapPlot);
        hold on;
        plot(t,x(1)+x(2)*t + x(3)*t.^2 + x(4)*t.^3, colors(pointsIteration));
		
% Save GA/path results
        numGenerations = Output.generations;
        fitnessValue = Fval;
        outputData(j,:) = [cSpaceID pointSetID x solutionLength dcoll maxJerk numGenerations fitnessValue gaLengthTime PopulationSize startPt endPt];
        end
		
% Save results file
    save(strcat(simdir,'/', fileName, '.txt'), 'outputData', '-ASCII', '-append');

    end
	
% Save path figure
    figname = strcat(simdir, '/' ,cSpaceFilenames(cSpaceIteration,:), '.fig');
    saveas(figure(cSpaceIteration), figname);
end

strcat('_Done',int2str(runNumber),'\n')

end