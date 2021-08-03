
% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = solution()
[consoleout, A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = evalc('solomiya_solution()');
end

function [A1, A2, A3, A4, A5, A6, A7, A8, A9, A10, A11, A12, A13, A14, A15, A16, A17, A18] = solomiya_solution()
% your solution code goes here
% assign the variables you are asked to save here
%% Homework 2 - Problem 1 "Shooting Method"

close all; clear all; clc

%defining the tolerance level
tolerance = 10^(-4);

%colors of eigenfunction
col = ['r', 'b', 'c', 'g', 'k', 'm'];

%defining the intial shooting value
x_starting = 1;

%defining the intial conditions
y0 = 1;

%defining the span of the domain
L = 4;
xspan = -L:0.1:L;

%initializing the column vectors
eigen_vec = zeros(5, 1);
i_modes = zeros(5, 1000);
y = zeros(81, 5);

%creating a for loop for the 5 eigenvecs
for modes = 1:5
    
    epsilon = x_starting;
    delta_time = 0.1;
    
    %convergence loop for x
    for i = 1:1000
        
        y_guess = sqrt(L^2 - epsilon);
        y_initial = [y0; y_guess];
        
        %defining the function
        shooting_first = @(t,y) [y(2); (t^2 - epsilon)*y(1)];
        
        %using ode45 to solve the function
        [t, y] = ode45(shooting_first, xspan, y_initial);
        
        %boundary solution equation
        bc = y(end, 2) + sqrt(L^2 - epsilon)*y(end, 1);
        
        if abs(bc) < tolerance
            eigen_vec(modes) = epsilon;
            y_extract(:, modes) = y(:, 1);
            break;
        end
        
        if (-1)^(modes) * bc > 0
            
            %check to see if epsilon needs to be higher or lower
            epsilon = epsilon - delta_time/2;
            delta_time = delta_time/2;
        else
            
            epsilon = epsilon + delta_time;
            
        end
    end
    
    x_starting = epsilon + 0.1;
    
    %normalizing the vectors
    norm = trapz(t, y_extract(:,modes).^2);
    %getting hw solutions
    y_extract(:,modes) = y_extract(:,modes)/sqrt(norm);
%     plot(t, y(:,1)/sqrt(norm), col(modes));hold on
    
end

A1 = abs(y_extract(:,1));
A2 = abs(y_extract(:,2));
A3 = abs(y_extract(:,3));
A4 = abs(y_extract(:,4));
A5 = abs(y_extract(:,5));
A6 = eigen_vec;

% AA = [A1 A2 A3 A4 A5]
%% Homework 2 - Problem 2 "Direct Method"
% clc; clear variables; close all;

%Setting parameters like del_t, K and L
L = 4;
x = -(L) : 0.1 : L;
del_time = 0.1;
K = 1;

%Creating a matrix of 79 rows 79 colums
N = 79;
F = zeros(N, N);

%Creating first and last row
two_thirds = 2 / 3;
F(1,1) = two_thirds + (del_time ^ 2) * K * (x(2) ^ 2);
F(1,2) = -(two_thirds);
F(N,N - 1) = -(two_thirds);
F(N,N) = two_thirds + (del_time ^ 2) * K * (x(N + 1) ^ 2);

%Creating a Loop for all the CDFs
for k = 2 : N - 1
    F(k, k - 1) = -1;
    F(k, k) = 2 + ((del_time ^ 2) * K * (x(k+1) ^ 2));
    F(k, k + 1) = -1;
end

%Solve for eigen values and eigen vectors
[y, epsilon] = eig(F);

%finding epsilon
epsilon = epsilon / (del_time ^ 2);

%Sorting the eigen values
eig_values_sort = sort(epsilon(epsilon > 0));
lowest5 = eig_values_sort(1 : 5);

%creating for loop to calculate outside boundary points
for i = 1:5
    
    %getting the vectors with the lowest eigens.
    [~, columnIndex] = find(epsilon == lowest5(i));
    running_column = y(:, columnIndex);
    
    x_1 = running_column(1);
    x_2 = running_column(2);
    
    %solving for inital and end boundary
    formula = (3 + 2 * del_time * sqrt(L^2 - lowest5(i)));
    y0 = (L * x_1 - x_2) / (formula);
    yn = (-(L) * running_column(N) + running_column(N - 1))/(formula);
    
    %adding the intial values and ending values to the eigen vector
    running_full_eigen_vectors = [y0; running_column; yn];
    
    %normalizing the vectors
    norm = trapz(x, running_full_eigen_vectors.^2);
    final_out(:, i) = abs(running_full_eigen_vectors/sqrt(norm));
    
    %Building a plot
%     figure(2),plot(x,final_out(:,i)); hold on;
    
end

%Solutions A7 - A12
A7 = final_out(:, 1);
A8 = final_out(:, 2);
A9 = final_out(:, 3);
A10 = final_out(:, 4);
A11 = final_out(:, 5);
A12 = lowest5;

% % AAA = [A7 A8 A9 A10 A11]

% Homework 2 - Problem 3 "Shooting Method"

% close all; clear all; clc
%defining the tolerance level
tolerance = 10^(-4);

%colors of eigenfunction
col = ['r', 'b', 'c', 'g', 'k', 'm'];

%defining the intial shooting value
x_start = 1;

%defining the intial conditions
AZ = 1;
% y0 = 0.1;
% K = 1;
gam_one = 0.05;

%defining the span of the domain
L = 2;
xspan = -L:0.1:L;

%initializing the column vectors
z_modes = zeros(2, 1);
y_extract = zeros(41, 2);

%creating a for loop for the 2 eigenvecs
for modes = 1:2
    
    %initialzing starting values
    ep_sil_on = x_start;
    delta_time = 1/55;
    
    %convergence loop for x
    for i = 1:1000
        %         epsilon = x_start;
        %         delta_time = 0.01;
        y_guess = sqrt(L^2 - ep_sil_on)*AZ;
        y_initial = [AZ; y_guess];
        
        
        %defining the function
        shooting_third = @(t,y) [y(2); ((gam_one*(y(1).^2) + t^2 - ep_sil_on)*y(1))];
        
        %using ode45 to solve the function
        [t, y] = ode45(shooting_third, xspan, y_initial);
        
        %computing the norm
        norm = trapz(t, y(:,1).^2);
        %          A = A/sqrt(norm);
        %boundary solution equation
        
        
        if abs(norm - 1) < tolerance
            break
        else
            AZ = AZ/sqrt(norm);
            
        end
        
        %boundary solution equation
        bc = y(end, 2) + sqrt(L^2 - ep_sil_on)*y(end, 1);
        
        if abs(bc) < tolerance
            
            z_modes(modes) = ep_sil_on;
            y_extract(:, modes) = abs(y(:, 1));
            break;
        end
        
        if (-1)^(modes + 1) * bc > 0
            ep_sil_on = ep_sil_on + delta_time;
            %check to see if epsilon needs to be higher or lower
            
        else
            
            ep_sil_on = ep_sil_on - delta_time/2;
            delta_time = delta_time/2;
            
        end
    end
    
    x_start = ep_sil_on + 0.1;
    
    %calculating norms again for eps
    norm = trapz(t, abs(y_extract(:,modes)).^2);
    
    %getting right solutions for hw
    y_extract(:,modes) = abs(y_extract(:,modes)/sqrt(norm));
    %      plot(t, y_extract(:,modes), col(modes));hold on
    
end


A13 = abs(y_extract(:,1));
A14 = abs(y_extract(:,2));
A15 = z_modes;

% AA = [A13 A14]
% Homework 2 - Problem 3 "Shooting Method"

% close all; clear all; clc
%defining the tolerance level
tolerance = 10^(-4);

%colors of eigenfunction
col = ['r', 'b', 'c', 'g', 'k', 'm'];

%defining the intial shooting value
x_start = 1;

%defining the intial conditions
AZ = 1; % defining A
% y0 = 0.1;
% K = 1;
gam_one = -0.05;

%defining the span of the domain
L = 2;

xspan = -L:0.1:L;

%initializing the column vectors
z_modes = zeros(2, 1);
y_extract = zeros(41, 2);

%creating a for loop for the 2 eigenvecs
for modes = 1:2
    
    %initialzing starting values
    ep_sil_on = x_start;
    delta_time = 1/11;
    
    %convergence loop for x
    for i = 1:1000
        %         epsilon = x_start;
        %         delta_time = 0.01;
        y_guess = sqrt(L^2 - ep_sil_on)*AZ;
        y_initial = [AZ; y_guess];
        
        
        %defining the function
        shooting_third = @(t,y) [y(2); ((gam_one*(y(1).^2) + t^2 - ep_sil_on)*y(1))];
        
        %using ode45 to solve the function
        [t, y] = ode45(shooting_third, xspan, y_initial);
        
        %computing the norm
        norm = trapz(t, y(:,1).^2);
        %          A = A/sqrt(norm);
        %boundary solution equation
        
        
        if abs(norm - 1) < tolerance
            break
        else
            AZ = AZ/sqrt(norm);
            
        end
        
        %boundary solution equation
        bc = y(end, 2) + sqrt(L^2 - ep_sil_on)*y(end, 1);
        
        if abs(bc) < tolerance
            
            z_modes(modes) = ep_sil_on;
            y_extract(:, modes) = abs(y(:, 1));
            break;
        end
        
        if (-1)^(modes + 1) * bc > 0
            
            ep_sil_on = ep_sil_on + delta_time;
            %check to see if epsilon needs to be higher or lower
            
        else
            
            ep_sil_on = ep_sil_on - delta_time/2;
            delta_time = delta_time/2;
            
        end
    end
    
    x_start = ep_sil_on + 0.1;
    
    %calculating norms again for eps
    norm = trapz(t, abs(y_extract(:,modes)).^2);
    
    %getting right solutions for hw
    y_extract(:,modes) = abs(y_extract(:,modes)/sqrt(norm));
    %      plot(t, y_extract(:,modes), col(modes));hold on
    
end


A16 = abs(y_extract(:,1));
A17 = abs(y_extract(:,2));
A18 = z_modes;

end
