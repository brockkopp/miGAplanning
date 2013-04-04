function plotWorld( world, armPos)
    
    [xDim yDim] = size(world);
    worldAxes = [0 xDim 0 yDim];

    figure(1); clf; hold on; title('My World')
        axis (worldAxes); 
        image(100*(1-world)');
        colormap(gray);

    plot([armPos(1,1) armPos(2,1)],[armPos(1,2) armPos(2,2)],'-b','LineWidth',5);
    plot([armPos(2,1) armPos(3,1)],[armPos(2,2) armPos(3,2)],'-b','LineWidth',5);

end

