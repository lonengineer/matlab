function [I, R] = romberg_integration(f, a, b, n, tol, options)
% romberg_integration: Implements Romberg integration to estimate the definite integral
% of a function f(x) over the interval [a, b].
%
% Inputs:
%   f       - Function handle for the function to be integrated.
%   a       - Lower limit of integration.
%   b       - Upper limit of integration.
%   n       - Maximum number of iterations (default: 100).
%   tol     - Desired tolerance for convergence (default: 1e-6).
%   options - Optional structure with additional parameters:
%               showDetailed: Logical flag to display detailed results (default: false).
%
% Outputs:
%   I - Final estimated value of the integral.
%   R - Romberg table showing intermediate results.
%
% Example:
%   % Estimate the integral of exp(-x^2) from 0 to 1 with detailed output
%   f = @(x) exp(-x.^2);
%   a = 0;
%   b = 1;
%   [I, R] = romberg_integration(f, a, b, [], 1e-8, showDetailed=true);
%   disp(I);
%   disp(R);
%
%   % Another example with custom iterations and default tolerance
%   [I, R] = romberg_integration(@(x) sin(x), 0, pi, 5);
%   fprintf('Integral estimate: %.10f\n', I);

arguments
    f
    a
    b
    n = 100
    tol = 1e-6
    options.showDetailed = false
end

maxIterations = n; % Maximum number of iterations
R = zeros(maxIterations, maxIterations); % Romberg table initialization
converged = false; % Flag to track convergence

% Step 1: Compute R(1,1) using the trapezoidal rule
R(1,1) = 0.5 * (b - a) * (f(a) + f(b));

for k = 2:maxIterations
    % Step 2: Compute the composite trapezoidal rule for 2^(k-1) subintervals
    num_subintervals = 2^(k-1); % Number of subintervals
    h = (b - a) / num_subintervals; % Subinterval width
    x = a + h * (1:2:num_subintervals-1); % Midpoints of subintervals
    R(k,1) = 0.5 * R(k-1,1) + h * sum(f(x));

    % Step 3: Richardson extrapolation
    for j = 2:k
        R(k,j) = R(k,j-1) + (R(k,j-1) - R(k-1,j-1)) / (4^(j-1) - 1);
    end

    % Check for convergence
    if abs(R(k,k) - R(k-1,k-1)) < tol
        converged = true;
        I = R(k,k);
        R = R(1:k, 1:k); % Trim the table to converged size
        break;
    end
end

% Handle non-convergence
if ~converged
    I = R(maxIterations, maxIterations);
    R = R(1:maxIterations, 1:maxIterations); % Return full table
    warning('Romberg integration did not converge to the desired tolerance within %d iterations.', maxIterations);
end

% Display detailed results if requested
if options.showDetailed
    fprintf('Estimated integral: %.8f\n', I);
    disp('Romberg table:');
    format long
    disp(R);
    format short
end
end
