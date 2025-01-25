function [df, E] = bd_diff(f, x0, h, options)
% bd_diff: Backward difference method for numerical differentiation
%
% Inputs:
%   f       - Symbolic function to differentiate
%   x0      - Point at which to evaluate the derivative
%   h       - Step size for backward difference
%   options - Optional name-value pair:
%               showDetailed: Logical flag to display detailed results (default: false)
%
% Outputs:
%   df - Backward difference approximation of f'(x) at x0
%   E  - Estimated truncation error bound
%
% Example:
%   syms x;
%   f = exp(x) + x^3;
%   [df, E] = bd_diff(f, 2, 0.001, 'showDetailed', true);
%   fprintf('Derivative estimate: %.5f\nError bound: %.2e\n', df, E);

arguments
    f
    x0
    h
    options.showDetailed (1,1) logical = false
end

% Compute function values at x0 and x0-h
F = subs(f, [x0, x0-h]);
df = (F(1) - F(2)) / h;

% Calculate second derivative and find maximum absolute value
f_dbl_prime = diff(f, 2);
t = x0-h:1e-5:x0;
sym_var = symvar(f);
f_dbl_prime_vals = double(subs(f_dbl_prime, sym_var, t));
[V, I] = max(abs(f_dbl_prime_vals));

% Calculate error bound
E = abs(V * h / 2);

% Display detailed results if requested
if options.showDetailed
    fprintf('\n--- Backward Difference Calculation Details ---\n');
    fprintf('Step size: \t\t%.5f\n', h);
    fprintf('x0: \t\t\t%.5f\n', x0);
    fprintf('f(x0): \t\t\t%.5f\n', double(F(1)));
    fprintf('f(x0-h): \t\t%.5f\n', double(F(2)));
    fprintf('\nDerivative estimate: \t%.8f\n', df);
    fprintf('Max |f''''(z)|: \t\t%.8f at z = %.5f\n', V, t(I));
    fprintf('Error bound: \t\t%.5e\n', E);
    fprintf('---------------------------------------------\n');
end
end
