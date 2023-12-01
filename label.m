function label(Title, Xlabel, Ylabel)
%LABEL is used to label a plot
% LABEL(Title, Xlabel, Ylabel), Title is title of plot,
% Xlabel is x-axis label,
% Ylabel is y-axis label
% LABEL with no arguments sets title as PLOT, x-axis as Time (s), y-axis as
% Amplitude

    arguments
        Title (1,:) {mustBeText} = 'Plot'
        Xlabel (1,:) {mustBeText} = 'Time (s)'
        Ylabel (1,:) {mustBeText} = 'Amplitude'
    end
    title(Title);
    xlabel(Xlabel);
    ylabel(Ylabel);
end
