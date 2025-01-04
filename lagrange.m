function L = lagrange(points, value, options)
%LAGRANGE Function
% This function computes the Lagrange polynomial for a given set of points,
% evaluates it at a specific value (optional), and provides additional
% features like error bound calculation and visualization.
%
% Syntax:
%   L = lagrange(points)
%   L = lagrange(points, value)
%   L = lagrange(points, value, options)
%
% Inputs:
%   points (:, 2) - A matrix where each row represents a data point (x, y).
%                   Example: [1, 2; 2, 4; 3, 6]
%
%   value (optional) - The x-coordinate at which to evaluate the polynomial.
%                      Default is '' (no evaluation). It simply return
%                      polynomial.
%
%   options (optional) - A struct with the following fields:
%       showPlot (logical) - Whether to display a plot of the polynomial
%                            and the points. Default: false.
%       xaxis (1, :)       - A vector defining the x-axis range for plotting.
%                            Default: auto-generated based on points.
%       calcError (logical) - Whether to calculate the interpolation error.
%                             Default: false.
%       expression (sym)   - The original symbolic function for error
%                            calculation. Required if calcError is true.
%
% Outputs:
%   L - The symbolic Lagrange polynomial or its evaluated value if 'value' 
%       is specified.
%
% Example Usage:
%   % Define points
%   points = [1, 2; 2, 4; 3, 6];
%   
%   % Interpolation with plot
%   L = lagrange(points, '', struct('showPlot', true));
%
%   % Interpolation with value evaluation
%   L_at_2_5 = lagrange(points, 2.5);
%
%   % Interpolation with error bound calculation
%   expr = @(x) x^2; % Original function
%   options = struct('calcError', true, 'expression', expr);
%   L = lagrange(points, calcError=true, expression=expr);

arguments
    points (:, 2)
    value = ''
    options.showPlot = false
    options.xaxis (1,:) = [points(1,1):(1/10)*(points(2,1)-points(1,1)): points(end,1)]
    options.calcError (1,1) logical = false
    options.expression 
    options.showDetailed (1,1) logical = false
end

    syms x;
    n = size(points, 1);  % Number of points based on rows of input
    L = 0;                % Initialize the Lagrange polynomial

    % Loop through each point to construct the Lagrange basis polynomials
    for i = 1:n
        % Initialize the basis polynomial for the current i
        L_i = 1;
        for j = 1:n
            if j ~= i
                % Multiply terms for the basis polynomial L_i
                L_i = L_i * (x - points(j, 1)) / (points(i, 1) - points(j, 1));
            end
        end
        % Add the term to the Lagrange polynomial
        L = L + points(i, 2) * L_i;
    end

    % Display the Lagrange polynomial for verification
    disp('Lagrange Polynomial, L(x):');
    L = simplify(L);
    L = vpa(L, 5);
    disp(L);

    if options.calcError
        pol = options.expression;
        e = 1;
        for ii = 1:n
            pol = diff(pol);
            e = e * (x - points(ii, 1));
        end
        
        range = [points(1,1), points(end,1)];
        [PMax, ~] = GMaxMin(pol, range);
        [eMax, ~] = GMaxMin(e, range);
        if options.showDetailed
            fprintf("%i derivative of original function:", n)
            pol
            fprintf("Polynomial equation for error:")
            e
            fprintf("Substituting %f in derivative and %f in error polynomial to maximize their values.\n", PMax(1), eMax(1))
        end
        EB = PMax(2)/factorial(n)*eMax(2);
        fprintf("Error bound is %f for the interval [%.2f, %.2f]", EB, points(1,1), points(end,1));
    end
    if options.showPlot
        x = options.xaxis;
        plot(points(:,1), points(:,2), 'o', LineWidth=3, MarkerSize=5, DisplayName="Points")
        hold on
        plot(x, subs(L,x), LineWidth=3, DisplayName="Interpollated")
        label 'Lagrange Polynomial Plot' 'x' 'L(x)'
        legend show; grid on;
        hold off
    end
    if ~ isempty(value), L = double(subs(L, value)); end
end
% 
% function L = lagrange(points, value, options)
% arguments
%     points (:, 2)
%     value = ''
%     options.showPlot = false
%     options.xaxis (1,2) = [-size(points, 1)*10, size(points, 1)*10]
% end
% 
%     syms x;
%     n = size(points, 1);  % Number of points based on rows of input
%     L = 0;                % Initialize the Lagrange polynomial
% 
%     % Loop through each point to construct the Lagrange basis polynomials
%     for i = 1:n
%         % Initialize the basis polynomial for the current i
%         L_i = 1;
%         for j = 1:n
%             if j ~= i
%                 % Multiply terms for the basis polynomial L_i
%                 L_i = L_i * (x - points(j, 1)) / (points(i, 1) - points(j, 1));
%             end
%         end
%         % Add the term to the Lagrange polynomial
%         L = L + points(i, 2) * L_i;
%     end
% 
%     % Display the Lagrange polynomial for verification
%     disp('Lagrange Polynomial, L(x):');
%     L = simplify(L);
%     L = vpa(L, 5);
%     disp(L);
% 
%     if options.showPlot
%         fplot(L, options.xaxis, Marker="o", DisplayName="Data Points", MarkerSize=8, Linewidth=1.5,Color='r');
%         hold on
%         fplot(L, options.xaxis, LineWidth=2.5, DisplayName="Lagrange Polynomial", Color='b');
%         label("Lagrange Polynomial Plot", "x", "L(x)")
%         grid on
%         hold off
%         legend show
%     end
%     if ~ isempty(value), L = double(subs(L, value)); end
% end
