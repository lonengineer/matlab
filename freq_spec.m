function [f w] = freq_spec(signal,sampling_rate,axis_limit)
%FREQ_SPEC plots frequency spectrum of a provided signal
    % FREQ_SPEC(signal,samping_rate),   signal is the signal whose spectrum is to be plotted and sampling_rate is number of samples taken
    % This function fails if the sampling rate is low and frequency is very high (comparitively)
    % This function plots spectrum with respect to w (rad/s) and not f
    % This function return vectors [f,w] to plot frequency spectrum of a signal
    % axis_limit limits the x axis to axis_limit from -axis_limit

f = abs(fftshift(fft(signal)/length(signal)));
w = 2*pi.*(linspace(-sampling_rate/2, sampling_rate/2,length(f)));
if(nargout==0)
    plot(w,f,LineWidth=2)
    title('Frequency Spectrum');
    ylabel('M(\omega)','FontWeight','bold')
    xlabel('\omega','FontWeight','bold')
    if nargin >= 3
        xlim([-axis_limit axis_limit]);
    end
end
end
