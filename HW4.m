% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5] = solution()
[consoleout, A1, A2, A3, A4, A5] = evalc('solomiya_solution()');
end

function [A1, A2, A3, A4, A5] = solomiya_solution()

% your solution code goes here

% clear all; close all; clc

%Problem 1 HW 4
% Creating P value in x & y directions
p = 64;

% total size of matrix
r = p * p;

%creating dxdy
dxdy = (10 - (-10))/p;

%creating a vector of zeros and ones
k0 = zeros(r, 1);
k1 = ones(r, 1);

%creating a copy of the ones and zeros vectors
k2 = k1;
k4 = k0;

%Creating matrix A
%overwriting t^th value with zeros and ones
for w = 1:p
    
    k2(p * w) = 0;
    k4(p * w) = 1;
    
end

%Shifting to correct the positions in k3 & k5
k3(2:r, 1) = k2(1:r-1, 1);
k3(1, 1) = k2(r, 1);
k5(2:r, 1) = k4(1:r-1, 1);
k5(1, 1) = k4(r, 1);

%Creating matrix A
% placing diagonal elements with 1/dxdy^2
A = spdiags([k1 k1 k5 k2 -4*k1 k3 k4 k1 k1], ...
    [-(r-p) -p -p+1 -1 0 1 p-1 p (r-p)],r,r);
A(1, 1) = 2;
A = A/(dxdy^2);

%Creating matrix B
% placing diagonal elements with 1/2*dxdy
B = (1/(2*dxdy)) * spdiags([k1 -k1 k1 -k1], ...
    [-(r-p) -p p (r-p)],r,r);

%Creating matrix C
% placing diagonal elements with 1/2*dxdy
C = (1/(2*dxdy)) * spdiags([k5 -k2 k3 -k4], ...
    [-p+1 -1 1 p-1],r,r);

%defining intial grid size
v = 0.001;

%spatial domain of x and y
D = 10 - (-10);

%accounting for periodicity
x_2 = linspace(-D/2, D/2, p+1);
x_1 = x_2(1:p);

%accounting for periodicity
y_2 = linspace(-D/2, D/2, p+1);

%y-vector
y_1 = y_2(1:p);

%setting up 2D intial conditions
[x_2, y_2] = meshgrid(x_1, y_1);

%defining t-span
d_t = 0.5;
t_span = 0 : d_t : p/16;

%generating a Gaussian matrix
z1 = exp(-x_2.^2 - (1/D)*(y_2.^2));

%reshapping into a vector
z = reshape(z1, [p^2 1]);

equation = @(t1, o)(-B*(A\o)).*(C*o) + (C*(A\o)).*(B*o) + v.*A*o;

%Solving for A\b
tic;

%solving the method
[t_1, y_1] = ode45(equation, t_span, z);

%recording the time to solve A\b equation
finish_1 = toc;

%solution for A1 A\b
A1 = y_1;

% LU method
[L, U, P] = lu(A);

%solving the method
equation_1 = @(t2, o)(-B*(U\(L\(P*o)))).*(C*o) + (C*(U\(L\(P*o)))).*(B*o) + v.*A*o;

tic;

[t_2, y_2] = ode45(equation_1, t_span, z);

%recording the time to solve the LU method
finish_2 = toc;

%solution for LU method
%LU method is much faster
A2 = y_2;

%Solving part c using FFT

[f_x,f_y] = meshgrid((2*pi/D)*[0:(p-1)/2 -p/2:-1],(2*pi/D)*[0:(p-1)/2 -p/2:-1]);
f_x(1) = 10e-6;
f_y(1) = f_x(1);

tic;

%solving the method
[t3, y_5] = ode45(@(t_3,o) equation3(v,o, A, B, C, p, f_x,f_y), t_span, z);

finish_3 = toc;

%solution for A5 using FFT
A5 = y_5;


%solving equation with bicgstab
[t_5, y_5] = ode45(@(t_5, o) equation5(v, o, A, B, C, z), t_span, z);
A3 = y_5;


%solving gmres equation with ode45
[t_4, y_4] = ode45(@(t_4, o) equation6(v, o, A, B, C, z), t_span, z);

%answer for A3
A4 = y_4;

%creating a function for FFT
    function fft = equation3(v, o, A, B, C, p, f_x,f_y)
        
        
        h = reshape(o, p, p);
        fft = -(B*(reshape(real(ifft2(-fft2(h)./(f_x.^2 + f_y.^2))), p^2, 1))).*(C*o) + (C*(reshape(real(ifft2(-fft2(h)./(f_x.^2 + f_y.^2))), p^2, 1))).*(B*o) + v*A*o;
    end


%gmres function command parameters
    function part2 = equation6(v, o, A, B, C, z)
        
        tolerance = 10^-5;
        [e, first, second, third] = gmres(A, o, [], tolerance, 4096);
        
        
        part2 =(-B*(e)).*(C*o) + (C*(e)).*(B*o) + v.*A*o;
        
    end

%bicgstab function command parameters
    function part = equation5(v, o, A, B, C, z)
        
        tolerance = 10^-5;
        [e1, first, second, third] = bicgstab(A, o, tolerance, 4096);
        
        
        part =(-B*(e1)).*(C*o) + (C*(e1)).*(B*o) + v.*A*o;
        
    end

end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)