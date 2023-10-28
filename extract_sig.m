function sig = extract_sig(sum_signal,sampling_rate,type,fc,order)
%EXTRACT_SIG extracts signal from a composit signal using butterworth filter, either high or low
% EXTRACT_SIG(sum_signal,sampling_rate,type,fc,order),  sum_signal is the composit signal, 
% sampling_rate is number of samples, 
% type is type of filter (high or low), 
% fc is the cutt off frequency of filter,
% order is order of filter used for extracting signal. 
% (default order is 4)

Wn = fc / (sampling_rate / 2);     % Normalized cutoff frequency
if nargin < 5
    order = 4;
end
[b, a] = butter(order,Wn , type);
sig = filter(b,a,sum_signal);
end