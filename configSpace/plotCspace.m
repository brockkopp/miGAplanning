function plotCspace( cSpace )

    [xDim yDim] = size(cSpace);
    cSpaceAxes = [0 xDim 0 yDim];
    
    figure(2); clf; hold on; title('Configuration Space')
        axis (cSpaceAxes); 
        image(100*(1-cSpace)');
        colormap(gray);
%     plot(armAng(1,1),armAng(1,2),'r+','MarkerSize',10);
%     plot(armAng(2,2),armAng(2,2),'b+','MarkerSize',10);
end

