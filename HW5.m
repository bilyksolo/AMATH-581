% Homework MATLAB template file
% Your main file should be named "solution.m" and it should be saved as UTF-8 file.

function [consoleout, A1, A2, A3] = solution()
 [consoleout, A1, A2, A3] = evalc('solomiya_solution()'); 
end

function [A1, A2, A3] = solomiya_solution()
    
%     clear all; close all;

    %setting up the intial parameters
    %for part a of homework problem
    L = 10 - (-10); %domain
    n = 64;

    %number of spirals
    m = 1; 

    %creating the span
    delta_t = 0.5;
    t_span = 0 : delta_t : n / 16;
    b = 1;

    %D - D1 and D2 from homework instructions
    D = 0.1;

    %accounting for periodicity
    x1_y1 = linspace(-L / 2, L / 2, n + 1);
    x1 = x1_y1(1 : n);
    y1 = x1;
    [X, Y] = meshgrid(x1, y1);

    %creating a grid
    formula = (2 * (pi/L)) * [0 : ((n/2 - 1)) ((-n/2)) : -1];
    formula(1) = 10^-6;

    [X1, Y1] = meshgrid(formula, formula);

    %multiplying the formula by a negative to get the right equation
    formula1 = (-1) * (X1.^2 + Y1.^2);

    %initializing u vector, v vector, uv vector
    uu = tanh(sqrt(X.^2 + Y.^2)).*cos(m * angle(X + 1i*Y) - (sqrt(X.^2 + Y.^2)));
    vv = tanh(sqrt(X.^2 + Y.^2)).*sin(m * angle(X + 1i*Y) - (sqrt(X.^2 + Y.^2)));

    %reshaping into a vector for u and v
    %creating a uv vector
    uv_final = [reshape(fft2(uu), n^2, 1); reshape(fft2(vv), n^2, 1)];

    %solving the equation
    [tt, yy] = ode45(@(tt2, yy2) periodic_fft(yy2, n, formula1, D, b), t_span, uv_final);

    %output for A1 and A2
    A1 = real(yy); %real solution
    A2 = imag(yy); %imaginary solution

    %solving part b
    %of homework 5 using n=30
    n = (n/2) - 1;

    [x, y] = cheb(n - 1);

    %rescaling
    x2 = (2/L).*x; y2 = (L/2).*y;

    x4 = x2 * x2;

    %creating the first and last row of zeros
    %for Dirichlet boundaries
    x4(1, :) = zeros(1, n);
    x4(n, :) = zeros(1, n);

    %creating the grid
    [X, Y] = meshgrid(y2, y2);

    %initializing u vector, v vector, uv vector
    %for part b of the problem
    uu1 = tanh(sqrt(X.^2 + Y.^2)).*cos(m * angle(X + 1i*Y) - (sqrt(X.^2 + Y.^2)));
    vv1 = tanh(sqrt(X.^2 + Y.^2)).*sin(m * angle(X + 1i*Y) - (sqrt(X.^2 + Y.^2)));

    %reshaping into a vector for u and v
    %creating a uv vector
    uv_final2 = [reshape(uu1, n^2, 1); reshape(vv1, n^2, 1)];

    %using kron command to create the eqation from the meshgrid
    equation = (kron(x4, eye(length(x4)))) + (kron(eye(length(x4)), x4));

    %solving the equation
    [tt1, yy2] = ode45(@(tt3, yy3) equation_chebyshev(yy3, n, equation, D, b), t_span, uv_final2);

    %output for A3
    A3 = yy2;

    %sample code from
    %class files
    %for Chebyshev equation
    function [DD1, x5] = cheb(H)
        if H == 0, DD1 = 0;
            x5 = 1;
            return;
        end
        x5 = cos(pi * (0:H)/H)';

        ccc = [2; ones(H-1, 1); 2].*(-1).^(0:H)';

        X = repmat(x5, 1, H + 1);

        d_X_1 = X - X';

        % off-diagonal entries
        DD1  = (ccc*(1./ccc)')./(d_X_1+(eye(H+1)));

        % diagonal entries
        DD1  = DD1 - diag(sum(DD1'));
    end

    %function provided from the 
    %homework assignment
    function partb = equation_chebyshev(vector, n, lap_lacian, D, b)

        %creating u and v vectors
        VVV = vector((n^2 + 1) : end);
        UUU = vector(1 : n^2);

        %inputting the vectors into a formula
        partb = [(1 - (UUU.^2 + VVV.^2)).*UUU - (-b*(UUU.^2 + VVV.^2)).*VVV + D*lap_lacian*UUU;...
            (-b*(UUU.^2 + VVV.^2)).*UUU + (1 - (UUU.^2 + VVV.^2)).*VVV + D*lap_lacian*VVV];
    end

    %Periodic boundary condition
    %function using FFT method
    function fft_periodic = periodic_fft(vector, n, lap_lacian, D, b)

        %reshaping the u and v vectors
        VVV1 = reshape(vector((n^2 + 1):end), n, n);
        UUU1 = reshape(vector(1:n^2), n, n);

        %using ift2 and real commands and then reshaping into vectors
        U_U = reshape((real(ifft2(UUU1))), n^2, 1);
        V_V = reshape((real(ifft2(VVV1))), n^2, 1);

        %inputting the vectors into a formula
        fft_periodic = [reshape((fft2(reshape(((1 - (U_U.^2 + V_V.^2)).*U_U - (-b*(U_U.^2 + V_V.^2)).*V_V), n, n))...
            + D*lap_lacian.*UUU1), n^2, 1); reshape((fft2(reshape(((-b*(U_U.^2 + V_V.^2)).*U_U...
            + (1 - (U_U.^2 + V_V.^2)).*V_V), n, n)) + D*lap_lacian.*VVV1), n^2, 1)];
    end
end

% your extra functions, if you need them, can be in other files (don't forget to upload them too!)