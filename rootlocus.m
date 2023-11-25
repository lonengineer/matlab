function rootlocus(G)
 % ROOTLOCUS  Evans root locus.
 % 
 %    rootlocus(SYS) computes and plots the root locus of the single-input,
 %    single-output LTI model SYS and visulasises how poles moves at different
 %    values of K. The root locus plot is used to analyze the negative feedback loop
 % 
 %                      +-----+
 %          ---->O----->| SYS |----+---->
 %              -|      +-----+    |
 %               |                 |
 %               |       +---+     |
 %               +-------| K |<----+
 %                       +---+
 % 
 %    and shows the trajectories of the closed-loop poles when the feedback 
 %    gain K varies from 0 to Inf.  rootlocus automatically generates a set of 
 %    positive gain values that produce a smooth plot.  

ni = nargin;
if ni==0,
    if no~=0,  error('Missing input argument(s).'),  end
    eval('exresp(''rlocus'');')
    return
end
[z,k]=rlocus(G);
label("Root Locus","Real Axis (seconds^{-1})","Imaginary Axis (seconds^{-1})")
hold on
num = asymtotes(G);
z=z(:,1:length(z)-(num+1));
z=z(:);
x_l = min(real(z));
x_h = max(real(z));
y_l = min(imag(z)); y_h = max(imag(z));

if (x_h == 0 && x_l == 0)
    x_h = 5; x_l = -5;
end
if (y_h == 0 && y_l == 0)
    y_h = 5; y_l = -5;
end

x_lim = [x_l+(x_l*0.05) x_h+(x_h*0.05)];
y_lim = [y_l+(y_l*0.05) y_h+(y_h*0.05)];
asymtotes(G);
xlim(x_lim)
ylim(y_lim)
for i = 1:length(k)-num
     if k(i)~=0 & k(i)~=Inf 
        sys = feedback(series(tf(k(i)),G),tf(1));
        [~, den] = tfdata(sys);
        r = roots(den{1});
        hold on
        plot(real(r),imag(r),'rx',LineWidth=1)
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
if nargout==0
    angleplot(theta,real(c),imag(c))
end
end

function angleplot(angle,x,y)
lineLength2 = 65535;
hold on;
for i = 1:numel(angle)
  e = x + lineLength2*cosd(angle(i));
  h = y + lineLength2*sind(angle(i));
  plot([x e], [y h],LineWidth=1,Color=[.5 .5 .5],LineStyle="--");
end
hold off;
end
