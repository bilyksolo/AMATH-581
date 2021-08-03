% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15] = solution()
 [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15] = evalc('solomiya_solution()'); 
end

function [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15] = solomiya_solution()
%     close all; clear all; clc
    
    %%  1. A1-A3  Euler  
    % Initial conditions that will be used 
    % for each âˆ†t = 2^âˆ’2, 2^âˆ’3, 2^âˆ’4, . . . , 2^âˆ’8
    allTs = [2^(-2), 2^(-3), 2^(-4), 2^(-5), 2^(-6), 2^(-7), 2^(-8)];
    size = length(allTs);
    A2 = zeros(1,size);
    y0 = pi / sqrt(2);
    
    eulerf = @(t,y) -3*y.*sin(t);
    
    for i = 1 : size % 7
        % Assign running t
        dt = allTs(i);
        
        % Creating t span
        t = 0:dt:5;
        
        % Initializing an empty array of t
        l = length(t);
        
        y = zeros(1, l);
        
        % Assigning initial condition
        y(1)= y0;
        
        % apply required formula to each y
        for n = 1 : l-1
            y(n+1)= y(n) + dt * eulerf(t(n),y(n));
        end
        
        % find the current iteration error value
        Exact=@(t,y)y0*exp(3*(cos(t)-1));
        Error = (abs(Exact(t)-y));
        
        % add this mean error value to the matrix A2
        err = mean(Error);
        A2(i) = err;
    end
    
    % Tha last error, using 2^âˆ’8
    %     
    %     A1 = A2(:,end);
    A1 = y';

    % The slop of the best fit line
    pol = polyfit(log(allTs), log(A2), 1);
    A3 = pol(1, 1);


    %%  1. A4-A6  Heun
    allTs = [2^(-2), 2^(-3), 2^(-4), 2^(-5), 2^(-6), 2^(-7), 2^(-8)];
    size = length(allTs);
    A5 = zeros(1,size);
    y0 = pi / sqrt(2);
    
    func = @(t,y) -3*y.*sin(t);
    
    for i = 1 : size % 7
        % Assign running t
        dt = allTs(i);
        
        % Creating t span
        t = 0:dt:5;
        
        % Initializing an empty array of t
        l = length(t);
        
        y = zeros(1, l);
        
        % Assigning initial condition
        y(1)= y0;
        
        % apply required formula to each y
        for n = 1 : l-1
            k1 = func(t(n), y(n));
            k2 = func(t(n) + dt, y(n) + dt * k1);
            
            y(n+1) = y(n) + (dt / 2) * (k1 + k2);
        end
        
        % find the current iteration error value
        Exact=@(t,y)y0*exp(3*(cos(t)-1));
        Error = (abs(Exact(t)-y));
        
        % add this mean error value to the matrix A2
        err = mean(Error);
        
        % for testing purposes we can show error
        % disp(err);
        
        A5(i) = err;
    end
    
    % Tha last error, using 2^âˆ’8
    A4 = y';
    %     A4 = A5(:,end);
    
    % The slop of the best fit line
    pol = polyfit(log(allTs), log(A5), 1);
    A6 = pol(1, 1);
    
    %%  2a. A7 Van Der Polâ€™s equation
    % Initial conditions that will be used
    allEs = [0.1 1.0 20];
    size = length(allEs);
    A7 = [];
    
    for i = 1 : size % 3
        % Assign running epsilon
        runningEpsilon = allEs(i);
        
        func = @(t,y) [y(2); runningEpsilon*((1-y(1)^2)*y(2))-y(1)];
        [t, y] = ode45(func, [0:0.5:32], [sqrt(3);1]);
        
        firstColumn = y(:,1);
        A7 = [A7, firstColumn];
    end
    
    %%  2b. A8, A9, A10 Van Der Polâ€™s equation
    allTs = [10^(-4) 10^(-5) 10^(-6) 10^(-7) 10^(-8) 10^(-9) 10^(-10)];
    tspan = [0,32];
    y0 = [2;(pi^2)];
    
    totalA8 = [];
    totalA9 = [];
    totalA10 = [];
    
    funct1 = @(t,y) [y(2); ((1-y(1)^2)*y(2))-y(1)];
     
    for i = allTs
        opt1 = odeset('AbsTol',i,'RelTol',i);
        
        % A8
        [t, y] = ode45(funct1, tspan, y0, opt1);
        
        firstColumn = t(:, 1);
        d = mean(diff(firstColumn));
        
        totalA8 = [totalA8, d];
    end
    
    for i = allTs
        opt1 = odeset('AbsTol',i,'RelTol',i);
        
        % A9
        [t, y] = ode23(funct1, tspan, y0, opt1);
        
        firstColumn = t(:, 1);
        d = mean(diff(firstColumn));
        
        totalA9 = [totalA9, d];
    end
    
    for i = allTs
        opt1 = odeset('AbsTol',i,'RelTol',i);
        
        % A10
        [t, y] = ode113(funct1, tspan, y0, opt1);
        
        firstColumn = t(:, 1);
        d = mean(diff(firstColumn));
        
        totalA10 = [totalA10, d];
    end
    
    % The slop of the best fit line
    pol = polyfit(log(totalA8), log(allTs), 1);
    A8 = pol(1, 1);
    
    % The slop of the best fit line
    pol = polyfit(log(totalA9), log(allTs), 1);
    A9 = pol(1, 1);
    
    % The slop of the best fit line
    pol = polyfit(log(totalA10), log(allTs), 1);
    A10 = pol(1, 1);
    
%     function dydt = funct1(t, y)
%         dydt = [y(2); ((1-y(1)^2)*y(2))-y(1)];
%     end

    %%  3. A11, A12, A13, A14, A15 Neuron
    % setting intial conditions
    y_initial =[0.1; 0.1; 0; 0];
    t = 0:0.5:100;
    AA = zeros(length(t),4,5);
    d12_d21 = [0 0; 0 0.2; -0.1 0.2; -0.3 0.2; -0.5 0.2];
    
    %creating for loop for the neurons function
    for j = 1:5
        [t, y] = ode15s(@(t, y) funcODE(y, d12_d21(j,:)), t, y_initial);
        AA(:,:,j) = y;
    end
    
    %answer output for A11 through A15
    A11 = AA(:,:,1);
    A12 = AA(:,:,2);
    A13 = AA(:,:,3);
    A14 = AA(:,:,4);
    A15 = AA(:,:,5);    

        function dy = funcODE(y, d)
            % setting parameters
            a1 = 0.05;
            a2 = 0.25;
            b = 0.01;
            c = 0.01;
            I = 0.1;

            dy = zeros(4,1);
            dy(1) = -y(1).^3 + (1+a1).*y(1).^2 - a1.*y(1) - y(3) + I + d(1).*y(2);
            dy(2) = -y(2).^3 + (1+a2).*y(2).^2 - a2.*y(2) - y(4) + I + d(2).*y(1);
            dy(3) = b.*y(1) - c.*y(3);
            dy(4) = b.*y(2) - c.*y(4);
            
        end
end

