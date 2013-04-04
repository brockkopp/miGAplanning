function [ armPos ] = updateArmPos( armLen, armAng )
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

armPos = zeros(2,2);
armAng(1) = degtorad(armAng(1));
armAng(2) = degtorad(armAng(2));

armPos(1,:) = [0+armLen(1)*cos(armAng(1)) 0+armLen(1)*sin(armAng(1))];
armPos(2,:) = [armPos(1,1)+armLen(2)*cos(armAng(2)) armPos(1,2)+armLen(2)*sin(armAng(2))];

end

