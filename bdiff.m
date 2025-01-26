function yint = bdiff(x, y, xx, options)
% bdiff: Backward difference interpolation polynomial constructor and evaluator
%
%   Constructs a Newton backward difference interpolation polynomial from 
%   equally spaced data points and evaluates it at specified locations.
%
% Inputs:
%   x       - Column vector of independent variable values (must be equally spaced)
%   y       - Column vector of dependent variable values (same length as x)
%   xx      - Optional evaluation points (scalar or vector)
%   options - Name-value arguments:
%       showTable  - Display divided difference table [false|true]
%       showPoly   - Display symbolic polynomial [false|true]
%       showPlot   - Plot data and interpolant [false|true]
%
% Outputs:
%   yint    - Interpolated values at xx (if xx provided), 
%             or symbolic polynomial (if xx empty)
%
% Example 1: Basic interpolation
%   x = [1; 2; 3; 4];
%   y = [2; 5; 10; 17];
%   p = bdiff(x, y);
%   fprintf('Interpolation polynomial:\n');
%   disp(p)
%
% Example 2: Full-featured usage
%   x = linspace(0, pi, 5)';
%   y = sin(x);
%   xx = pi/4;
%   yi = bdiff(x, y, xx, 'showTable', true, 'showPoly', true, 'showPlot', true);
%   fprintf('sin(π/4) ≈ %.8f\nError: %.2e\n', yi, abs(yi - sin(pi/4)));
%
% Notes:
%   - Requires equally spaced x-values
%   - Polynomial construction uses symbolic math toolbox
%   - First divided difference is Δf = f(x_i) - f(x_{i-1})

    arguments
        x (:,1)  % Independent variable (x)
        y (:,1)  % Dependent variable (y)
        xx = []   % Points to interpolate (optional)
        options.showTable (1,1) logical = false  % Display the difference table
        options.showPoly (1,1) logical = false   % Display the interpolating polynomial
        options.showPlot (1,1) logical = false   % Plot the function and interpolant
    end

    % Ensure that x and y have the same length
    n = length(x); temp = x;
    if length(y) ~= n
        error('x and y must be the same length');
    end

    % Initialize the backward divided difference table
    b = zeros(n, n);
    b(:,1) = y(:); % Assign dependent variable to the first column of b

    % Compute the backward divided differences
    for j = 2:n
        for i = n:-1:j
            b(i, j) = (b(i, j-1) - b(i-1, j-1)); % Backward difference
        end
    end

    % Show the difference table if option is enabled
    if options.showTable
        headers = ["x", "Δ f"];
        for j = 1:n-1
            headers = [headers, sprintf("Δ^%i f", j+1)];
        end

        T = array2table([x b], 'VariableNames', string(headers));
        disp('Backward Difference Table:');
        disp(T);
    end

    % Use the finite divided differences to construct the interpolation polynomial
    h = x(2) - x(1); % Step size
    xn = x(end); % Last x value (for backward formula)
    syms s;  % Use symbolic variable for constructing the polynomial
    poly = b(n,1); % Start with f(xn)
    
    % Add terms using backward differences from last row of table
    for j = 1:n-1
        % Get backward difference coefficient from last row
        coeff = b(n, j+1); 
        % Calculate binomial coefficient with s + j - 1 choose j
        binomial = nchoosek(s + j - 1, j);
        % Add term to polynomial
        poly = poly + coeff * binomial;
    end
    poly_s = poly;
    % Substitute s with (x - xn)/h and expand
    syms x
    poly_x = subs(poly, s, (x - xn)/h);
    poly_x = expand(poly_x);
    
    % Display the polynomial if option is enabled
    if options.showPoly
        disp('Interpolating Polynomial P(s):');
        disp(poly_s);
        disp('Interpolating Polynomial P(x):');
        disp(poly_x);
    end
    
    % Expand the polynomial
    poly_x = expand(poly_x);
    
    % Evaluate the polynomial at xx (if provided) or return the polynomial
    if ~isempty(xx)
        yint = double(subs(poly_x, x, xx));
    else
        yint = poly_x;
    end

    % Plot the function and the interpolating polynomial if option is enabled
    if options.showPlot
        x = temp;
        plot(x, y, 'o', 'LineWidth', 2, 'DisplayName', 'Data Points');
        hold on;
        x_vals = linspace(min(x), max(x), 1000);
        plot(x_vals, double(subs(poly_x, x_vals)), 'LineWidth', 3, 'DisplayName', 'Interpolated Polynomial');
        legend show;
        grid on;
        xlabel('x');
        ylabel('P(x)');
        title('Backward Difference Interpolation');
    end
end
