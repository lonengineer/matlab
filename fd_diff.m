function [df, E] = fd_diff(f, x0, h, options)
% fd_diff: Forward difference derivative approximation with error estimation
%
% Inputs:
%   f       - Symbolic function to differentiate
%   x0      - Point at which to compute derivative
%   h       - Step size for forward difference
%   options - Optional name-value pair:
%               showDetailed: Logical flag to display detailed results (default: false)
%
% Outputs:
%   df - Forward difference estimate of f'(x0)
%   E  - Error bound estimate using maximum of |f''(z)| in [x0, x0+h]
%
% Example:
%   syms x;
%   f = sin(x) + x^2;
%   [df, E] = fd_diff(f, 0.5, 0.001, 'showDetailed', true);
%   fprintf('Estimate: %.5f\nError bound: %.2e\n', df, E);

arguments
    f
    x0
    h
    options.showDetailed (1,1) logical = false
end

% Compute function values at x0 and x0+h
F = subs(f, [x0, x0+h]);
df = (F(2) - F(1))/h;

% Calculate second derivative and find maximum
f_dbl_prime = diff(f, 2);
t = x0:1e-5:x0+h;
sym_var = symvar(f);
f_dbl_prime_vals = double(subs(f_dbl_prime, sym_var, t));
[V, I] = max(abs(f_dbl_prime_vals));

% Calculate error bound
E = abs(V*h/2);

% Detailed output formatting
if options.showDetailed
    fprintf('\n--- Forward Difference Calculation Details ---\n');
    fprintf('Step size: \t\t%.5f\n', h);
    fprintf('x0: \t\t\t%.5f\n', x0);
    fprintf('f(x0): \t\t\t%.5f\n', double(F(1)));
    fprintf('f(x0+h): \t\t%.5f\n', double(F(2)));
    fprintf('\nDerivative estimate: \t%.8f\n', df);
    fprintf('Max |f''''(z)|: \t\t%.8f at z = %.5f\n', V, t(I));
    fprintf('Error bound: \t\t%.5e\n', E);
    fprintf('---------------------------------------------\n');
end
end
