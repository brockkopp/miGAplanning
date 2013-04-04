function [ funval ] = AKfitness( x, S, E, obsWeight, lengthWeightFactor, lineResolution )
HIGH_PUNISHMENT = 10000;

%
startPt = S;
endPt = E;
A = x(1);
B = x(2);
C = x(3);
D = x(4);
E = x(5);

%Import grid from workspace
obstacleGrid = evalin('base','obsGrid');
[xDim yDim] = size(obstacleGrid);

MAX_Y = yDim;
MIN_Y = 0;

%% TO DO: GET RID OF SYMBOLIC INTEGRATION

%y = sqrt(1 + B^2);
% intvalue=int(y,t,startPt(1),endPt(1));
% funval=double(intvalue);


% syms t; 
% t = startPt(1):0.1:endPt(1);
% y = sqrt((B^2 + 1) + (4*B*C)*t + (4*C^2)*t*t);

y2 = @(t) sqrt(1 + (B + 2*C*t + 3*D*t.^2 + 4*E*t.^3).^2); % (B^2 + 1) + (4*B*C)*t + (4*C^2)*t.^2);
dist = integral(y2, startPt(1), endPt(1));

% Punish for Obstacle
obstaclePunishment = 0;
for i=startPt(1):lineResolution:endPt(1)
    y = A + B*i + C*i^2 + D*i^3 + E*i^4;
    if (y > MAX_Y || y < MIN_Y)
        %obstaclePunishment = HIGH_PUNISHMENT;
    obstaclePunishment = obstaclePunishment + obsWeight;
    elseif (obstacleGrid(i,ceil(y)) == 1) %if within obstacle
        obstaclePunishment = obstaclePunishment + obsWeight;
    end
end

funval = dist*lengthWeightFactor + obstaclePunishment;
end
