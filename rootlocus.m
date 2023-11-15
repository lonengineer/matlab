function rootlocus(G)
% rlocus(G)
% fig = gcf;
% x_lim=fig.CurrentAxes.XLim;
% y_lim=fig.CurrentAxes.YLim;
% close
[~,k]=rlocus(G);
label("Root Locus","Real Axis (seconds^{-1})","Imaginary Axis (seconds^{-1})")
num = asymtotes(G);
% xlim(x_lim)
% ylim(y_lim)
hold on
for i = 1:length(k)-num
     if k(i)~=0 & k(i)~=Inf 
        sys = feedback(series(tf(k(i)),G),tf(1));
        [~, den] = tfdata(sys);
        r = roots(den{1});
        hold on
        plot(real(r),imag(r),'rx',LineWidth=2)
        pause(1/50)
     end
end
hold off
end

function n = asymtotes(sys)
[num,denum]=tfdata(sys);
p = roots(denum{1});
z = roots(num{1});
n = length(p)-length(z);
centroid = @(zero,pole)(sum(pole)-sum(zero))/(length(pole)-length(zero));
phi=@(l)(180+360*(l-1))/(length(p)-length(z));
c=centroid(z,p);
theta=[];
for i = 1:(length(p)-length(z))
    temp=phi(i);
    theta = [theta temp];
end
angleplot(theta,real(c),imag(c))
end

function angleplot(angle,x,y)
lineLength2 = 30;
hold on;
for i = 1:numel(angle)
  e = x + lineLength2*cosd(angle(i));
  h = y + lineLength2*sind(angle(i));
  plot([x e], [y h],LineWidth=1,Color=[.5 .5 .5],LineStyle="--");
end
hold off;
end
