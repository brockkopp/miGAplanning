function [startPt, endPt] = generatePoints(cSpace)

[xLim, yLim] = size(cSpace);

invalid = true;
while(invalid)
    startPtx = ceil((xLim/2)*rand());
    startPty = ceil((yLim)*rand());
    
    invalid = (cSpace(startPtx, startPty) == 1);
end

invalid = true;
while(invalid)
    endPtx = ceil((xLim/2) + (xLim/2)*rand());
    endPty = ceil((yLim)*rand());
    
    invalid = (cSpace(endPtx, endPty) == 1);
end

startPt = [startPtx startPty];
endPt = [endPtx endPty];

end
