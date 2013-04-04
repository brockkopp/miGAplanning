% Integrate exp(x)*sin(3*x) from x=2.0 to 8.7
% Define x as a symbol
syms x
% Assigning the function to be differentiated
y=exp(x)*sin(3*x);
% Assigning the lower limit
a=2.0;
% Assigning the upper limit
b=8.7;

%% DISPLAYING INPUTS

disp('INPUTS')
func=['  The function is to be integrated is ' char(y)];
disp(func)
fprintf('  Lower limit of integration, a= %g',a)
fprintf('\n  Upper limit of integration, b= %g',b)
disp(' ')
disp(' ')

%% THE CODE

% Finding the integral using the int command
% Argument 1 is the function to be integrated
% Argument 2 is the variable with respect to which the
%    function is to be integrated – the dummy variable
% Argument 3 is the lower limit of integration
% Argument 4 is the upper imit of integration
intvalue=int(y,x,a,b);
intvalue=double(intvalue);

%% DISPLAYING OUTPUTS

disp('OUTPUTS')
fprintf('  Value of integral is = %g',intvalue)
disp(' ')