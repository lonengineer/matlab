function Wp = pcf(sys)
syms wp real
[n d] = numden(sys);
z = roots(sym2poly(n));
z = abs(z);
p = roots(sym2poly(d));
p = abs(p);
eq = 0;
for i = 1:length(z)
    if z(i)==0
        eq = eq+90;
    else
        eq = eq+atand(wp/z(i));
    end
end
for i = 1:length(p)
    if p(i)==0
        eq = eq-90;
    else
        eq = eq-atand(wp/p(i));
    end
end
Wp = solve(eq==-180,wp)
Wp = double(Wp);
end