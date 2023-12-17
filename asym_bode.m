function  [y x y_phase x_phase] = asym_bode(wn,n,lim)
arguments
    wn, n = 1, lim = [10e-3 10e5];
end
x{length(wn)}=0; y{length(wn)}=0;
b(length(wn))=0;
m = n*20;
for i = 1:length(wn)
    if wn(i)==0
        k=1; 
    else
        k = wn(i);
    end
    b(i) = 0-m*log10(k);

    x{i} = logspace(log10(lim(1)),log10(lim(2)),1000);
    y{i} = m*log10(x{i})+b(i);
    if wn(i)~=0
        y{i}(x{i}<=k) = 0;
    end
    
end
subplot(211);
for i = 1:length(wn)
    semilogx(x{i},y{i},LineWidth=1);
hold on
% set(gca,'xtick',[])
end
title 'Asymtotic Bode Plot'
ylabel 'Gain (dB)'

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
m_phase = n*45;
x_phase{length(wn)}=0; y_phase{length(wn)}=0;
for i = 1:length(wn)
    x_phase{i} = logspace(log10(lim(1)),log10(lim(2)),1000);
    if wn(i)==0
        y_phase{i} = ones(size(x_phase{i}))*2*m_phase;
    else
        k = log10(wn(i));
        k = [10^(k-1) 10^(k+1)];
        b(i) = 0-m_phase*log10(k(1));
        y_phase{i} = m_phase*log10(x_phase{i})+b(i);
        y_phase{i}(x{i}<k(1)) = 0;
        y_phase{i}(x{i}>k(2)) = 2*m_phase;

    end
    
    
end
subplot(212);

for i = 1:length(wn)
    semilogx(x_phase{i},y_phase{i},LineWidth=1);
hold on
end
% set(gca,'xtickMode', 'auto')
ylabel Phase
xlabel 'Frequency \omega_n'
hold off
end
