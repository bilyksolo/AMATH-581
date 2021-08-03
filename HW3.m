% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3, A4, A5] = solution()
[consoleout, A1, A2, A3, A4, A5] = evalc('solomiya_solution()');
end

function [A1, A2, A3, A4, A5] = solomiya_solution()

% keep these lines to load Fmat and permvec
% DO NOT submit the mat files to Gradescope
load Fmat.mat
load permvec.mat

% your solution code goes here

% assign the variables you are asked to save here
%     A1 = 0;
%     A2 = 0;
%     A3 = 0;
%     A4 = 0;
%     A5 = 0;

%% Problem #2 of HW3
%getting the length of permvec file
%the order of permve: 7 11 3 14 4 16 5 15 2 6 1 10 13 9 8 12
%needs to be        : 1 2 3 4 5 6  7 8  9 10 11 12 13 14 15 16
span_of_per_mvec = length(permvec);

%initializing the center of Fmat
center_of_Fmat = Fmat(161:240, 161:240);

%size of the cube
dimentions = 80 / 4;

%creating a matrix of zeros 80 by 80
Fmat_zeros = zeros(20 * 4);

%creating row and column of orignal permvec and
%the updated version
[r_o_w_20, c_o_l_20] = ind2sub([16 / 4, 16 / 4], 1 : span_of_per_mvec);
[r_o_w_10, c_o_l_10] = ind2sub([16 / 4, 16 / 4], permvec);
E = (dimentions - 1);

%creating a for loop to look up orignal blocks
%and shifting them to an appropriate place
for w= 1:16
    
    %creating row and column with original dimentions
    A = r_o_w_20(w)*dimentions;
    B = c_o_l_20(w)*dimentions;
    
    %modifing the row and column
    C = r_o_w_10(w)*dimentions;
    D = c_o_l_10(w)*dimentions;
    
    %updating the matrix with correct cubes of picture
    Fmat_zeros(((A - E) : A), ((B - E) : B)) = center_of_Fmat(((C - E) : C), ((D - E) : D));
    
end

%initializing the beginning image
first_copy = Fmat;

%getting the first image
first_copy(161:240, 161:240) = Fmat_zeros;

%final image product
megan_and_harr = abs(ifft2(ifftshift(first_copy)));

%answers for the assignment
A4 = abs(first_copy);
A5 = megan_and_harr;

%% Problem 1 of HW3
% clear all; close all;   % clear all variables and figures

% code from page 55 from the notes
% Creating P value in x & y directions
p = 8;

% total size of matrix
t = p * p;

%creating dxdy
dxdy = (10 - (-10))/8;

%creating a vector of zeros and ones
k0 = zeros(t, 1);
k1 = ones(t, 1);

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
k3(2:t, 1) = k2(1:t-1, 1);

k3(1, 1) = k2(t, 1);

k5(2:t, 1) = k4(1:t-1, 1);

k5(1, 1) = k4(t, 1);

% placing diagonal elements with 1/dxdy^2
A_matrix_final = (1/(dxdy^2)) * spdiags([k1 k1 k5 k2 -4*k1 k3 k4 k1 k1], ...
    [-(t-p) -p -p+1 -1 0 1 p-1 p (t-p)],t,t);


% answer for A matrix
A1 = full(A_matrix_final);

% viewing the matrix structure using the spy command
% spy(A_matrix_final);

%Creating matrix B
% placing diagonal elements with 1/2*dxdy
B_matrix_final = (1/(2*dxdy)) * spdiags([k1 -k1 k1 -k1], ...
    [-(t-p) -p p (t-p)],t,t);


% answer for B matrix
A2 = full(B_matrix_final);

% viewing the matrix structure using the spy command
% spy(B_matrix_final);


%Creating matrix C
% placing diagonal elements with 1/2*dxdy

C_matrix_final = (1/(2*dxdy)) * spdiags([k5 -k2 k3 -k4], ...
    [-p+1 -1 1 p-1],t,t);

% Answer for C matrix
A3 = full(C_matrix_final);

% viewing the matrix structure using the spy command
% spy(C_matrix_final);


end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)