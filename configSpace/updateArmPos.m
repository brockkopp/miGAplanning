function [ armPos ] = updateArmPos( armBase, armLen, armAng )

armPos = zeros(3,2);
armAng(1) = degtorad(armAng(1));
armAng(2) = degtorad(armAng(2));

armPos(1,:) = armBase(:);
armPos(2,:) = [armPos(1,1)+armLen(1)*cos(armAng(1)) armPos(1,2)+armLen(1)*sin(armAng(1))];
armPos(3,:) = [armPos(2,1)+armLen(2)*cos(armAng(2)) armPos(2,2)+armLen(2)*sin(armAng(2))];

end

