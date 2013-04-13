%outputData(j,:) = 
%     [    1- cSpaceID,
%          2- pointSet
%          3- x - A
%          4- x - B
%          5- x - C
%          6- x - D
%          7- x - E
%          8- solutionLength 
%          9- numCollisions 
%         10- maxJerk 
%         11- numGenerations 
%         12- fitnessValue 
%         13- gaLengthTime 
%         14- PopulationSize 
%         15- startPt - X
%         16- startPt - Y
%         17- endPt - X
%         18- endPt - Y];

function [ output_args ] = getGoldStd( gaData, cSpace1, cSpace2, cSpace3, cSpace4 )

    cSpace = [cSpace1 cSpace2 cSpace3 cSpace4];

    points = gaData(:, 15:18);
    cSpaceId = gaData(:, 1);
    [numRows, ~] = size(points);

    pathLengths = ones(length(points), 1);

    if mod(numRows,numPts) ~= 0
        error('Length mismatch');
    else
        numTests = numRows / numPts;
    end

    results = zeros(numTests,5);

    for i=1:numTests-1
        
%% Load results
        set = zeros(5,length(gaData(1,:)));

        startIdx = i*numPts +1;
        set(1,:) = gaData(startIdx,:);

        endY = set(1,18);

        for j=2:4
            idx = startIdx + j;
            if(gaData(idx,18) ~= endY)
               error('Coordinates do note match') 
            else
                set(j,:) = gaData(startIdx,:);
            end
        end
        
        startPt = set(1,:)
        
        wfPathLength = wavefront(startPt, endPt, cSpace(set);
        
%%WAVEFRONT

    

    %           startPt = [points(i, 1) points(i, 2)];
    %           endPt   = [points(i, 3) points(i, 4)];

              startPt = [50 200];
              endPt = [250 250];
              pathLengths(i) = 
          end
        
%%Print Results
        results(i+1,1) = mean(set(:,8));  % Length Mean
        results(i+1,2) = mean(set(:,8));  % Length Std Dev
        results(i+1,3) = min(set(:,8));   % Length min
        results(i+1,4) = max(set(:,9));   % Max collisions
        results(i+1,5) = max(set(:,12));  % Max fitness

    end

    
%% Wavefront
    for i=1:numRows

    end
    
    


    save('wavefrontData.txt', 'pathLengths', '-ASCII');
    disp '_Done'
end

