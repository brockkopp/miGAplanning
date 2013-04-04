function [ cSpace, F ] = buildCspace( armBase, armLen, world, limits, F )

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
            if(alpha > 360)
                cSpace(alpha,beta) = cSpace(alpha-360,beta);
            else
                armPos = updateArmPos(armBase, armLen,[alpha beta]);

    %             if(mod(beta,45) == 0)
    %                 plotAll( world, armPos, cSpace )
    %                 F(frm) = getframe(figure(3));
    %                 frm = frm +1;
    %             end

                cSpace(alpha,beta) = 0;

                if(armPos(2,1) <= 1 || armPos(2,1) >= xDim-1 || ...
                   armPos(2,2) <= 1 || armPos(2,2) >= yDim-1 || ...
                   armPos(3,1) <= 1 || armPos(3,1) >= xDim-1 || ...
                   armPos(3,2) <= 1 || armPos(3,2) >= yDim-1)
                    cSpace(alpha,beta) = 1;
                else
                    for i=1:armLen(1)
                        x1 = round(armPos(1,1)+i*cos(degtorad(alpha))+1);
                        y1 = round(armPos(1,2)+i*sin(degtorad(alpha))+1);
                        if(world(x1,y1) == 1)
                            cSpace(alpha,beta) = 1;
                            break;
                        end
                    end
                    if(cSpace(alpha,beta) == 0)
                        for i=1:armLen(2)
                            x2 = round(armPos(2,1)+i*cos(degtorad(beta))+1);
                            y2 = round(armPos(2,2)+i*sin(degtorad(beta))+1);
                            if(world(x2,y2) == 1)
                                cSpace(alpha,beta) = 1;
                                break;
                            end
                        end
                    end %End arm2
                end %End Limits
            end %End if 360
        end %End For
    end %End alpha loop
    
    plotAll( world, armPos, cSpace );
    F(frm) = getframe(figure(3));
    close(h);
end

