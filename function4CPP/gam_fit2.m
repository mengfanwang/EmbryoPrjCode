rr = gamrnd(2,0.2, 10000,1);
p = gamfit(rr)


a0 = 0.5*(log(mean(rr)) - mean(log(rr)));
a = a0;
for i = 1:5
    d = (psi(a+0.0001) - psi(a-0.0001)) / 0.0002;
    a = 1 / (1/a + (mean(log(rr)) - log(mean(rr)) + log(a) - psi(a)) ...
        / (a*a*(1/a - d)));
end

b = mean(rr) / a;
disp([a b])