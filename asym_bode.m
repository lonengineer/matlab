function  [y x] = asym_bode(wn,n,lim)
arguments
    wn, n = 1, lim = [10e-3 10e9];
end
m = n*20;
for i = 1:length(wn)
    if wn(i)==0
        k=1; 
    else
        k = wn(i);
    end
    b(i) = 0-m*log10(k);

    x{i} = logspace(log10(lim(1)),log10(lim(2)));
    y{i} = m*log10(x{i})+b(i);
    if wn(i)~=0
        y{i}(x{i}<=k) = 0;
    end
    
end

for i = 1:length(wn)
    semilogx(x{i},y{i},LineWidth=2);
hold on
end
grid on
label 'Asymtotic Bode Plot' 'Frequency \omega_n' 'Gain (dB)'
hold off
end
