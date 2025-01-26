function L = lagrange(points, value, options)
%LAGRANGE Construct and evaluate Lagrange interpolation polynomial
%
%   Computes the Lagrange interpolating polynomial for given data points,
%   evaluates it at specified locations, and provides error estimation
%   with visualization options.
%
% Syntax:
%   L = lagrange(points)
%   L = lagrange(points, value)
%   L = lagrange(points, value, options)
%
% Inputs:
%   points  - (N×2 numeric) Matrix of [x, y] data points
%   value   - (numeric) Optional evaluation point(s)
%   options - Name-value arguments:
%       showPlot     - Plot polynomial and points [false|true]
%       xaxis        - Custom x-range for plotting [auto-generated]
%       calcError    - Calculate error bound [false|true]
%       expression   - (symbolic) Original function for error calculation
%       showDetailed - Display computation details [false|true]
%
% Outputs:
%   L - (symbolic/double) Polynomial or evaluated value
%
% Examples:
%   % Basic interpolation
%   pts = [1 1; 2 4; 3 9];
%   P = lagrange(pts);
%   
%   % Evaluate at specific point with error estimation
%   syms x;
%   L = lagrange(pts, 2.5, 'showDetailed', true, ...
%                'calcError', true, 'expression', x^2);
%
%   % Plot interpolation
%   lagrange(pts, [], 'showPlot', true);

arguments
    points (:,2) double
    value {mustBeNumeric} = []
    options.showPlot (1,1) logical = false
    options.xaxis (1,:) double = linspace(min(points(:,1)), max(points(:,1)), 1000)
    options.calcError (1,1) logical = false
    options.expression sym = sym([])
    options.showDetailed (1,1) logical = false
end

% Symbolic variable setup
syms x;
n = size(points, 1);
L_sym = sym(0);

% Construct Lagrange polynomial
for i = 1:n
    xi = points(i, 1);
    yi = points(i, 2);
    basis = sym(1);
    
    for j = 1:n
        if j ~= i
            xj = points(j, 1);
            basis = basis * (x - xj) / (xi - xj);
        end
    end
    L_sym = L_sym + yi * basis;
end

% Simplify and prepare output
L = simplify(L_sym);
if ~isempty(value)
    L = double(subs(L, x, value));
end

% Error bound calculation
if options.calcError
    if isempty(options.expression)
        error('Original function expression required for error calculation');
    end
    
    % Calculate nth derivative
    f_deriv = options.expression;
    for k = 1:n
        f_deriv = diff(f_deriv);
    end
    
    % Error term components
    error_term = sym(1);
    for k = 1:n
        error_term = error_term * (x - points(k,1));
    end
    
    % Manual maxima calculation
    x_vals = linspace(min(points(:,1)), max(points(:,1)), 1000);
    
    % Evaluate derivatives
    deriv_vals = double(subs(abs(f_deriv), x, x_vals));
    [max_deriv, idx_deriv] = max(deriv_vals);
    
    % Evaluate error terms
    error_vals = double(subs(abs(error_term), x, x_vals));
    [max_error, idx_error] = max(error_vals);
    
    % Calculate error bound
    error_bound = (max_deriv/factorial(n)) * max_error;
    
    % Detailed output
    if options.showDetailed
        fprintf('\n--- Error Calculation Details ---\n');
        fprintf('Maximum %dth derivative: %.4f at x = %.4f\n', ...
                n, max_deriv, x_vals(idx_deriv));
        fprintf('Maximum error term: %.12f at x = %.4f\n', ...
                max_error, x_vals(idx_error));
        fprintf('Error bound: %.4e\n', error_bound);
        fprintf('---------------------------------\n');
    end
end

% Visualization
if options.showPlot
    figure;
    plot(points(:,1), points(:,2), 'o', 'MarkerSize', 8, ...
        'LineWidth', 2, 'DisplayName', 'Data Points');
    hold on;
    y_vals = double(subs(L_sym, x, options.xaxis));
    plot(options.xaxis, y_vals, 'LineWidth', 2, ...
         'DisplayName', 'Interpolation');
    title('Lagrange Interpolation');
    xlabel('x');
    ylabel('P(x)');
    legend show;
    grid on;
    hold off;
end

% Display polynomial if no output requested
if nargout == 0
    disp('Lagrange Interpolating Polynomial:');
    (L_sym)
end
end
