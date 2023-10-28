function label(Title, Xlabel, Ylabel)
%LABEL is used to label a plot
% LABEL(Title, Xlabel, Ylabel), Title is title of plot,
% Xlabel is x-axis label,
% Ylabel is y-axis label
% LABEL with no arguments sets title as PLOT, x-axis as Time (s), y-axis as
% Amplitude
    if nargin<1
        Title = 'PLOT';
        Xlabel = 'Time (s)';
        Ylabel = 'Amplitude';
    elseif nargin<2
        Xlabel = 'Time (s)';
        Ylabel = 'Amplitude';
    elseif nargin<3
        Ylabel = 'Amplitude';
    end
    title(Title);
    xlabel(Xlabel);
    ylabel(Ylabel);
end