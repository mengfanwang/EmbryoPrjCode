function [mu, sigma] = ordStatApprox3sec(fg, bg, order_vec, alpha, beta)
% approximate the ~normal distribution using the fg and bg vector
% assuming there is only three parts exist for [0 1]
%
M = length(fg);
N = length(bg);
n = M+N;
a1 = n/M;
a2 = -n/N;
% |--alpha--|--beta--|--1-alpha-beta--|
beta = alpha + beta; % get the real second cutting point
%% mu: we use the original formular to calculate
J = zeros(1,n);
J(order_vec(:,1) <0) = -n/N;
J(order_vec(:,1) >0) = n/M;
u = 1/(n+1):1/(n+1):n/(n+1);
Fmu = J.*norminv(u);
mu = 1/2*(1/(n+1))*(2*sum(Fmu) - Fmu(1) - Fmu(end));
%% mu: we use the estimation: this will be lower
% mu = a2*(-normpdf(norminv(alpha))) + a1*normpdf(norminv(beta));
%% sigma we use the 3-sec estimation wit inf
%     inf_val = 1e3;
%     c1 = 2*a2^2*norminv(alpha)+2*a1*a2*(inf_val-norminv(beta));
%     t1 = f1(norminv(alpha));
%     c2 = 2*a1^2*inf_val;
%     t2 = f1(inf_val) - f1(norminv(beta));
%     c3 = -2*a2*a2;
%     t3 = f2(norminv(alpha));
%     c4 = -2*a1*a1;
%     t4 = f2(inf_val) - f2(norminv(beta));
% 
%     A = c1*t1+c2*t2+c3*t3+c4*t4;
%     B = (a2*f1(norminv(alpha)) + a1*(inf_val-f1(norminv(beta))))^2;
%% sigma we use the 3-sec estimation wit inf
c1 = 2*a2^2*norminv(alpha)+2*a1*a2*(-norminv(beta));
t1 = f1(norminv(alpha));
c3 = -2*a2*a2;
t3 = f2(norminv(alpha));
c4 = -2*a1*a1;
t4 = -0.5-f2(norminv(beta));

A = c1*t1+c3*t3+c4*t4;
B = (a2*f1(norminv(alpha)) + a1*(-f1(norminv(beta))))^2;%1.91

sigma = sqrt(A-B)/sqrt(n);

end

function y = f1(x)
y = x*normcdf(x)+normpdf(x);
end

function y = f2(x)
y=0.5*(normcdf(x)*x^2 - normcdf(x) + normpdf(x)*x);
end
