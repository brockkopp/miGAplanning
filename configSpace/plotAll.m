function [fig] = plotAll( world, armPos, cSpace )

    [WxDim WyDim] = size(world);
    worldAxes = [0 WxDim 0 WyDim];
    
    [CxDim CyDim] = size(cSpace);
    cSpaceAxes = [0 CxDim 0 CyDim];

    fig = figure(3); 
    clf;
        subplot(1,2,1);
            hold on;
            image(100*(1-world)');
            colormap(gray);
            plot([armPos(1,1) armPos(2,1)],[armPos(1,1) armPos(2,2)],'-b','LineWidth',5);
            plot([armPos(2,1) armPos(3,1)],[armPos(2,2) armPos(3,2)],'-r','LineWidth',5);
            axis (worldAxes); 
            title('Workspace','FontSize',14);
            xlabel('X-Position (cm)','FontSize',12)
            ylabel('Y-Position (cm)','FontSize',12)
        
        subplot(1,2,2);
            image(100*(1-cSpace)');
            colormap(gray);
            axis (cSpaceAxes);
            title('Configuration Space','FontSize',14);
            xlabel('Joint 1 Angle (deg)','FontSize',12,'Color','b')
            ylabel('Joint 2 Angle (deg)','FontSize',12,'Color','r')
    
end

