function [ fit ] = splineFitness( dna, cSpace, armAng, X )
    fit = 0;
    
    [~,yDim] = size(cSpace);
    offset = armAng(1,2) - (dna(4)*armAng(1,1)^4 + dna(3)*armAng(1,1)^3 + dna(2)*armAng(1,1)^2 + dna(1)*armAng(1,1));
    
    adj = offset*ones(size(X));
    
    Y = dna(4)*X.^4 + dna(3)*X.^3 + dna(2)*X.^2 + dna(1)*X + adj; 

    for x = armAng(1,1):armAng(2,1)
        if(Y(x) <= 0 || Y(x) >= yDim || cSpace(round(X(x)),round(Y(x))) == 1)
            fit = fit + 10;
        end
    end
    
    endDistance = abs(armAng(2,2) - Y(armAng(2,1)));
    fit = fit + endDistance;

end

