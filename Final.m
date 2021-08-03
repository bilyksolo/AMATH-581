% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4] = solution()
[consoleout, A1, A2, A3, A4] = evalc('solomiya_solution()');
end

function [A1, A2, A3, A4] = solomiya_solution()

clear all; close all;

%Final Project
%Bose-Einstein Condensation in 3D
n = 4 * 4;
L = 2 * pi;
delta_t = 0.5;
tspan = 0:delta_t:n/4;

%creating the grid space
w = linspace(-L/2, L/2, n+1);
X = w(1:n);
Y = w(1:n);
Z = w(1:n);
[xx, yy, zz] = meshgrid(X, Y, Z);

%reshaping x into a vector
x = reshape(xx, n^3, 1);

%reshaping y into a vector
y = reshape(yy, n^3, 1);

%reshaping z into a vector
z = reshape(zz, n^3, 1);

%defining the A and B matrices
AA = [-1 -1 -1];
BB = [1 1 1];

%creating the forier grid space
f = (2*pi/L)*[0:(n/2-1) (-n/2):-1];

%setting the tolerance
f(1) = 10^-6;

[DX, DY, DZ] = meshgrid(f, f, f);

%creating the laplacian
laplacian = (-1) * (DX.^2 + DY.^2 + DZ.^2);

%creating the intial condition for the equation
%intial condition with cose
psi_initial = cos(xx).*cos(yy).*cos(zz);

%creating a vector of intial condition
psi_vector = reshape(fftn(psi_initial), n^3, 1);

%solving the ode equation
[tt, yyy] = ode45(@(t, y1) rhs(t, y1, laplacian, AA, BB, x, y, z, n), tspan, psi_vector);

%solutions for A1 and A2
%A1 is for real part
A1 = real(yyy);

%A2 is for the imaginary part
A2 = imag(yyy);

% solving part b
% updating intial condition to sine
psi_initial_partb = sin(xx).*sin(yy).*sin(zz);

psi_vector_partb = reshape(fftn(psi_initial_partb), n^3, 1);

%solving the equation
[tt1, yy1] = ode45(@(t, y1) rhs(t, y1, laplacian, AA, BB, x, y, z, n), tspan, psi_vector_partb);

%A3 is for real part
A3 = real(yy1);

%A4 is for the imaginary part
A4 = imag(yy1);

    function bose = rhs(t, y1, laplacian, AA, BB, x, y, z, n)
        
        %creating psi
        psi1 = reshape(y1, n, n, n);
        
        %reshaping psi into a vector
        psi = reshape(ifftn(psi1), n^3, 1);
        
        %calculating the last part of the equation
        first = (AA(1).*sin(x).^2 + BB(1));
        second = (AA(2).*sin(y).^2 + BB(2));
        third = (AA(3).*sin(z).^2 + BB(3));
        part3 = (first.*second.*third);
        
        %second part of the equation
        part2 = (-conj(psi).*psi + part3).*psi;
        part = fftn(reshape(part2, n, n, n));
        
        part1 = (laplacian.*psi1)/2;
        
        %putting everything together and then
        %dividing by -1i to get the right solution
        equation = (part1 + part)/(-1i);
        bose = reshape(equation, n^3, 1);
        
    end
end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)