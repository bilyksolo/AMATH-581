% AMATH 581 Fall 2020
% Homework Gradescope Submission Practice
% Solomiya Bilyk

function [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = solution() 
	[consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = evalc('solomiya_solution()'); 
end 

function [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13] = solomiya_solution()    
	% Exercise 1 Building a Matrix
    A = [34 45; 17 6];
    A1 = A;

    % Exercise 2 Matrix Operations
    % Define variables
    A = [1 2;-1 1];
    B = [2 0; 0 2];
    C = [2 0 -3;0 0 -1];
    D = [1 2; 2 3; -1 0];
    x = [1 0].';
    y = [0 1].';
    z = [1 2 -1].';

    % Assign variables A2-A10
    A2 = A + B;
    A3 = 3 * x - 4 * y;
    A4 = A * x;
    A5 = B * (x - y);
    A6 = D * x;
    A7 = D * y + z;
    A8 = A * B;
    A9 = B * C;
    A10 = C * D;
    
    % Exercise 3 Root Finding
    
    % A11 is the column vector of x-values in the 
    % Newton method starting with the initial guess x(1) = -3
	A11 = []; 
         
    % Initialize known variables
    f = @(x) -x - cos(x);
    fd = @(x) -1 + sin(x);
    x0 = -3;    % lower
    x1 = 1;     % upper 
    tol = 10e-7; 
    
    x = x0;
    y = f(x);
    A11 = [A11; x'];
    
    while abs(y) > tol
        x = x - y / fd(x);
        y = f(x);
        A11 = [A11; x'];
    end
         
	% A12 is the column vector of midpoint (xmid) values 
    % in the bisection method for successive iterations
    A12 = []; 
    
    x_lower = x0;
    x_upper = x1;
    
    % find the mid value
    x_mid= (x_lower + x_upper) / 2.0;
    y = f(x_mid);
    
    % add beginning value to the matrix
    A12 = [A12, x_mid,];
    
    while abs(y) > tol        
        % check if we update x_lower of x_upper variable
        if ((f(x_mid) * f(x_upper)) < 0)
            x_lower = x_mid;
        else
            x_upper = x_mid;
        end
        
        % update the running mid value
        x_mid = (x_lower + x_upper) / 2.0;
        
        % reassing y to check if we succes tollerance in the next while
        % loop run
        y = f(x_mid);
        
        A12 = [A12; x_mid,];
    end
   
	% A13 is a 1x2 vector with the number of iterations 
    % for the Newton and bisection respectively as the two components    
    A13 = [length(A11), length(A12)];
end