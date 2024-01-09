function lldesign(sys,options)
arguments
    sys
    options.z, options.wn, options.leadz
    options.lagp, options.OSp, options.ts
end
syms s
damping = @(p_os)sqrt((log(p_os/100)^2)/(pi^2+log(p_os/100)^2));
nat_freq = @(zeta, settling_time) (4.6/(zeta*settling_time));

zeta =damping(options.OSp)
wn = nat_freq(zeta,options.ts)
desiredLocation = roots ([1 2*zeta*wn wn^2])

point = [real(desiredLocation(1)) abs(imag(desiredLocation(1)))]

[num den] = numden(sys);
psys = roots(sym2poly(den))
zsys = roots(sym2poly(num))

theta = zeros(1,length(psys)+length(zsys));
for i = 1:length(psys)
    if (real(psys(i))<point(1))
        x = abs(real(psys(i)))-abs(point(1));
        theta(i) = -(atand(point(2)/(x)));
    else
        x = abs(point(1))-abs(real(psys(i)));
        theta(i) = -(180-atand(point(2)/(x)));
    end
end

for i = 1:length(zsys)
    j = length(psys)+i;
    if (real(zsys(i))<point(1))
        x = abs(real(zsys(i)))-abs(point(1));
        theta(j) = atand(point(2)/(x));
    else
        x = abs(point(1))-abs(real(zsys(i)));
        theta(j) = 180-atand(point(2)/(x));
    end
end

leadz = options.leadz
if (leadz<point(1))
        x = abs(leadz)-abs(point(1));
        tz = atand(point(2)/(x));
    else
        x = abs(point(1))-abs(leadz);
        tz = 180-atand(point(2)/(x));
end

theta

% tz - tp = -180
tp = 180+sum(double(theta))+tz

leadp = -(point(2)/tand(tp))-abs(point(1))

K_lead = (s-leadz)/(s-leadp);
K_lead = vpa(K_lead)

end

