function [ cSpace, F ] = buildCspace( armLen, world, limits, F )

    h = waitbar(0,'Initializing');

    cSpace = 0.5*ones(limits(2)-limits(1)+1,limits(4)-limits(3)+1);
    [xDim yDim] = size(world);
    
    frm = 1;
%     fig = figure(3);
%     winsize = get(figure(3),'Position');
%     winsize(1:2) = [0 0];

    for alpha = limits(1):limits(2)
        waitbar(alpha/limits(2),h,'Processing...');
        
        for beta = limits(3):limits(4)
            
            armPos = updateArmPos(armLen,[alpha beta]);

%             if(mod(beta,45) == 0)
%                 plotAll( world, armPos, cSpace )
%                 F(frm) = getframe(figure(3));
%                 frm = frm +1;
%             end
            
            cSpace(alpha,beta) = 0;
            
            if(armPos(1,1) <= 0 || armPos(1,1) >= xDim || ...
               armPos(1,2) <= 0 || armPos(1,2) >= yDim || ...
               armPos(2,1) <= 0 || armPos(2,1) >= xDim || ...
               armPos(2,2) <= 0 || armPos(2,2) >= yDim)
                cSpace(alpha,beta) = 1;
            else
                for i=0:armLen(1)
                    x1 = round(i*cos(degtorad(alpha))+1);
                    y1 = round(i*sin(degtorad(alpha))+1);
                    if(world(x1,y1) == 1)
                        cSpace(alpha,beta) = 1;
                        break;
                    end
                end
                if(cSpace(alpha,beta) == 0)
                    for i=0:armLen(2)
                        x2 = round(armPos(1,1)+i*cos(degtorad(beta))+1);
                        y2 = round(armPos(1,2)+i*sin(degtorad(beta))+1);
                        if(world(x2,y2) == 1)
                            cSpace(alpha,beta) = 1;
                            break;
                        end
                    end
                end
            end
        end
    end
    plotAll( world, armPos, cSpace );
    F(frm) = getframe(figure(3));
    close(h);
end

