function yint = fdiff(x, y, xx, options)
% fdiff: Forward difference interpolation polynomial constructor and evaluator
%
%   Constructs a Newton forward difference interpolation polynomial from 
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
% Example 1: Basic polynomial interpolation
%   x = [1; 2; 3; 4];
%   y = [1; 4; 9; 16];
%   p = fdiff(x, y);
%   fprintf('Interpolation polynomial:\n');
%   disp(p)
%
% Example 2: Trigonometric function with full diagnostics
%   x = linspace(0, pi/2, 4)';
%   y = sin(x);
%   xx = pi/4;
%   yi = fdiff(x, y, xx, 'showTable', true, 'showPoly', true, 'showPlot', true);
%   fprintf('sin(π/4) ≈ %.8f\nError: %.2e\n', yi, abs(yi - sin(pi/4)));
%
% Notes:
%   - Requires equally spaced x-values
%   - Polynomial construction uses symbolic math toolbox
%   - First divided difference is Δf = f(x_{i+1}) - f(x_i)
%   - Follows Newton forward formula: 
%     P(x) = f(x₀) + sΔf + s(s-1)/2! Δ²f + ... where s = (x - x₀)/h

arguments
    x (:,1)
    y (:,1)
    xx = []
    options.showTable (1,1) logical = false
    options.showPoly (1,1) logical = false
    options.showPlot (1,1) logical = false
end

    % Compute the finite divided differences in the form of a
    % difference table
    n = length(x);
    if length(y) ~= n
        error('x and y must be the same length');
    end

    % Initialize the divided difference table
    b = zeros(n, n);
    % Assign dependent variables to the first column of b
    b(:,1) = y(:); 

    % Calculate the divided differences
    for j = 2:n
        for i = 1:n-j+1
            b(i, j) = (b(i+1, j-1) - b(i, j-1));
        end
    end

    if options.showTable
        headers = ["x", "Δ f"];
        for j = 1:n-1
            headers = [headers, sprintf("Δ^%i f",j+1)];
        end

        T = array2table([x b], 'VariableNames', string(headers));
        disp('Forward Difference Table:')
        disp(T)
    end

    % Use the finite divided differences to interpolate and construct the polynomial
    h = x(2)-x(1);
    temp = x; % Store the original x values
    syms s  % Use symbolic x for constructing the polynomial
    xt = 1;
    yint = b(1,1);
    poly = yint; % Start building polynomial with first term

    for j = 1:n-1
        xt = nchoosek(s,j);
        term = b(1, j+1) * xt;
        poly = poly + term; % Add each term to polynomial
    end
    poly_s = simplify(poly);
    syms x
    poly_x = subs(poly_s, s=(x-temp(1))/h);
    % Display the polynomial if option is enabled
    if options.showPoly
        disp('Interpolating Polynomial P(s):')
        poly_s
        disp('Interpolating Polynomial P(x):')
        poly_x
    end
    poly_x = expand(poly_x);
    if ~isempty(xx)
        yint = double(subs(poly_x, x, xx));
    else
        yint = poly_x;
    end

    if options.showPlot
        x = temp;
        plot(x, y,'o', LineWidth=2, DisplayName="Points")
        hold on
        x = x(1):1/100 * (x(2)-x(1)):x(end);
        plot(x, subs(poly_x, x), LineWidth=3, DisplayName="Interpollated")
        legend show; label 'Forward Difference Interpolation' 'x' P(x); grid on
    end

end
