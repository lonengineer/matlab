function K_lag = lag_design(sys,options)
arguments
    sys
    options.wn
    options.error = 0.1
    options.lagp
end
if ~isempty(options.wn)
    wn=options.wn
end
syms s
syms k z real

sys_ol = k*sys;

sys_cl = sys_ol/(1+sys_ol);
sys_cl = simplify(sys_cl);

[num den] = numden(sys_cl);

Kp = vpa(solve(subs(den,0)==wn^2,k));

[~, dens] = numden(sys);
temp = sym2poly(dens)

z = vpa(solve(2*z*wn==temp(end-1)))
subplot(223)
rlocus(symtotf(sys_ol))

E = 1/(1+sys_ol);
gain = solve(subs(E,s=0)==options.error,k)
subplot(221)
step(symtotf(sys_cl,Kp))
hold on
step(symtotf(E,Kp))
error = dcgain(symtotf(E,Kp))

lagp = options.lagp;

lag_gain = gain/Kp
lagz = lag_gain*lagp;
K_lag = (s-lagz)/(s-lagp)
sys_kol = k*((s-lagz)/(s-lagp))*sys;
sys_kcl = sys_kol/(1+(sys_kol));
subplot(222)
step(symtotf(sys_kcl,Kp))
hold on
E_k = 1/(1+sys_ol*K_lag);
step(symtotf(E_k,Kp))
[~, temp] = numden(sys_kcl)
subplot(224)
rlocus(symtotf(sys_kol/k))
dp = roots(sym2poly(subs(temp,k=Kp)))
hold on
plot(real(dp),imag(dp),'rx','MarkerSize',10,LineWidth=1)

;
end