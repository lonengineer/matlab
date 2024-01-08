function cotroller(sys,options)
damping = @(p_os)sqrt((log(p_os/100)^2)/(pi^2+log(p_os/100)^2));
nat_freq = @(zeta, settling_time) (4.6/(zeta*settling_time));

zeta =damping(4.5)
wn = nat_freq(zeta,9)
desiredLocation = roots ([1 2*zeta*wn wn^2])
point = [real(desiredLocation) abs(imag(desiredLocation))]

[num den] = numden(sys);
psys = roots(sym2poly(den))
zsys = roots(sym2poly(num))

theta = zeros(1,length(psys)+length(zsys));
for i = 1:length(psys)
    if (real(psys(i))<point(1))
        theta(i) = -(atand(point(2)/(abs(real(psys(i)))-abs(point(1)))))
    else
        theta(i) = -(180-atand(point(2)/(abs(point(1))-abs(real(psys(i))))))
    end
end

for i = 1:length(zsys)
    j = length(psys)+i;
    if (real(zsys(i))<point(1))
        theta(j) = atand(point(2)/(abs(real(zsys(i)))-abs(point(1))))
    else
        theta(j) = 180-atand(point(2)/(abs(point(1))-abs(real(zsys(i)))))
    end
end

syms tz tp real
% theta = vpa([theta tz -tp])
req = -180-sum(theta)-tz+tp
vpa(req)
end
