results = importdata('gaData.txt');

cSpace1 = importdata(strcat('../configSpace/',  'cSpace2.mat'));
cSpace2 = importdata(strcat('../configSpace/',  'cSpace3.mat'));
cSpace3 = importdata(strcat('../configSpace/',  'cSpace4.mat'));
cSpace4 = importdata(strcat('../configSpace/',  'obsGrid.mat'));

points = results(:, 14:17);
cSpaceId = results(:, 1);
[numRows numCols] = size(points);

pathLengths = ones(length(points), 1);
        
for i=1:numRows
      GoOn = true;
      if(cSpaceId(i) == 1)
          cSpace = cSpace1;
      elseif(cSpaceId(i) == 2)
          cSpace = cSpace2;
      elseif(cSpaceId(i) == 3)
          cSpace = cSpace3;
      elseif(cSpaceId(i) == 4)
          cSpace = cSpace4;
      else
          GoOn = false;
      end
      
      if(GoOn)
%           startPt = [points(i, 1) points(i, 2)];
%           endPt   = [points(i, 3) points(i, 4)];

          startPt = [50 200];
          endPt = [250 250];
          pathLengths(i) = wavefront(startPt, endPt, cSpace);
      end
end


save('wavefrontData.txt', 'pathLengths', '-ASCII');
disp '_Done'