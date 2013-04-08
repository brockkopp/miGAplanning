%outputData(j,:) = 
%     [    1- cSpaceID,
%          2- pointSet
%          3- x - A
%          4- x - B
%          5- x - C
%          6- x - D
%          7- solutionLength 
%          8- numCollisions 
%          9- maxJerk 
%         10- numGenerations 
%         11- fitnessValue 
%         12- gaLengthTime 
%         13- PopulationSize 
%         14- startPt - X
%         15- startPt - Y
%         16- endPt - X
%         17- endPt - Y];

clear all;

load('gaData_cub_600.mat');
load('gaData_cubJerk_600.mat');
load('wfData.mat');

numPts = 40;

gaData1 = gaData_cub_600;
gaData2 = gaData_cubJerk_600;
wfData1 = wfData;

if(size(gaData1) ~= size(gaData2))
    error('Size mismatch');
end


if mod(length(gaData1(:,1)),numPts) ~= 0
    error('Length mismatch');
else
    numTests = length(gaData1(:,1)) / numPts;
end

resLen = zeros(numTests,6);

diffMatrix = zeros(length(gaData1(:,1)),4);
shorterCnt = 0;
shorterCnt_noJerk = 0;
noColl = 0;
noColl_noJerk = 0;

for testNum=1:numTests

    set1 = zeros(5,length(gaData1(1,:)));
    set2 = zeros(5,length(gaData2(1,:)));
   
    startIdx = (testNum-1)*numPts +1;
    set1(1,:) = gaData1(startIdx,:);
    set2(1,:) = gaData2(startIdx,:);

    endY = set1(1,17);
    
    wfIdx = 5*(set1(1,1)-1) + set1(1,2);

    for j=2:(numPts -1)
        idx = startIdx + j;
        if(gaData1(idx,17) ~= endY || gaData2(idx,17) ~= endY)
           error('Coordinates do note match') 
        else
            set1(j,:) = gaData1(startIdx+j-1,:);
            set2(j,:) = gaData2(startIdx+j-1,:);
            
            diffMatrix(startIdx+j-1,1) = gaData1(startIdx+j-1,7) / wfData(wfIdx,3);
            diffMatrix(startIdx+j-1,2) = gaData2(startIdx+j-1,7) / wfData(wfIdx,3);
            
            if gaData2(startIdx+j-1,9) == 0
                diffMatrix(startIdx+j-1,3) = 0;
            else
                diffMatrix(startIdx+j-1,3) = gaData1(startIdx+j-1,9) / gaData2(startIdx+j-1,9);
            end
            
            
            if(gaData1(startIdx+j-1,7) < wfData(wfIdx,3))
                shorterCnt = shorterCnt +1;
            end
            if(gaData2(startIdx+j-1,7) < wfData(wfIdx,3))
                shorterCnt_noJerk = shorterCnt_noJerk +1;
            end
            
            if(gaData1(startIdx+j-1,8) == 0)
                noColl = noColl + 1;
            end
            if(gaData2(startIdx+j-1,8) == 0)
                noColl_noJerk = noColl_noJerk + 1;
            end
            
        end
    end

    resLen(testNum, 1) = mean(set1(:,7));  % Length Mean
    resLen(testNum, 2) = mean(set2(:,7));  % Length Mean
    resLen(testNum, 3) = wfData(wfIdx,3);  % Length Mean
    
    resLen(testNum, 5) = mean(set1(:,7))/wfData(wfIdx,3);  % Length Mean
    resLen(testNum, 6) = mean(set2(:,7))/wfData(wfIdx,3);  % Length Mean
    
    
    
%     resLen(testNum, 4) = std(set1(:,7));  % Length Std Dev
%     resLen(testNum, 5) = std(set2(:,7));  % Length Std Dev
%     resLen(testNum, 6) = 0;  % Length Std Dev
    
    results(testNum, 1) = max(set1(:,8));   % Max collisions
    results(testNum, 2) = max(set2(:,8));   % Max collisions
    results(testNum, 3) = max(set1(:,11));  % Max fitness
    results(testNum, 4) = max(set2(:,11));  % Max fitness  
end

% stdDiff = mean(diffMatrix(:,1))
% stdDiff_noJerk = mean(diffMatrix(:,2))
fprintf('\n\tRESULTS \n\n');
fprintf('Std:\t shorter  \t\t\t\t%d of %d\n', shorterCnt, length(gaData1(:,1)))
fprintf('Std:\t Mean difference:  \t\t%d\n', mean(diffMatrix(:,1)))
fprintf('Std:\t Col''n free:   \t\t\t%d of %d\n', noColl, length(gaData1(:,1)))
fprintf('Std:\t Mean max jerk comp:  \t%d\n', mean(diffMatrix(:,3)))
fprintf('Std:\t Mean comp time:  \t\t%d\n', mean(gaData1(:,12)))
fprintf('Std:\t Mean num generations:  %d\n', mean(gaData1(:,10)))
fprintf('Std:\t Std Dev generations:  %d\n', std(gaData1(:,10)))


fprintf('\n');
fprintf('NoJerk:\t shorter  \t\t\t\t%d of %d\n',shorterCnt_noJerk, length(gaData1(:,1)))
fprintf('NoJerk:\t Col''n free:   \t\t\t%d of %d\n', noColl_noJerk, length(gaData1(:,1)))
fprintf('NoJerk:\t Mean difference:  \t\t%d\n', mean(diffMatrix(:,2)))
fprintf('NoJerk:\t Mean comp time:  \t\t%d\n', mean(gaData2(:,12)))
fprintf('NoJerk:\t Mean num generations:  %d\n', mean(gaData2(:,10)))
fprintf('NoJerk:\t Std Dev generations:  %d\n', std(gaData2(:,10)))


