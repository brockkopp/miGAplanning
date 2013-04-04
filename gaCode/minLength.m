function [ funval ] = minLength( x, S, E )
%
startPt = S;
endPt = E;
A = x(1);
B = x(2);
C = x(3);
D = x(4);
E = x(5);

%% TO DO: GET RID OF SYMBOLIC INTEGRATION

%y = sqrt(1 + B^2);
% intvalue=int(y,t,startPt(1),endPt(1));
% funval=double(intvalue);


% syms t; 
% t = startPt(1):0.1:endPt(1);
% y = sqrt((B^2 + 1) + (4*B*C)*t + (4*C^2)*t*t);
% y = @(t) sqrt((B^2 + 1) + (4*B*C)*t + (4*C^2)*t.^2);
y = @(t) sqrt(1 + (B + 2*C*t + 3*D*t.^2 + 4*E*t.^3).^2); % (B^2 + 1) + (4*B*C)*t + (4*C^2)*t.^2);
funval = integral(y, startPt(1), endPt(1));
end
