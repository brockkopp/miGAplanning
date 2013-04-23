obsGrid = cSpace    
mapPlot = figure; 
    hold on;
    [xDim yDim] = size(obsGrid);
    axis([0 xDim 0 yDim]);
    image(100*(1-obsGrid)');
    colormap(gray);