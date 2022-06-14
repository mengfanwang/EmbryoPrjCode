function [mu, sigma] = ordStatApproxKsec(fg, bg)
% approximate the ~normal distribution using the fg and bg vector
% assuming there is only K parts exist for [0 1]
%
if isempty(fg) || isempty(bg)
    mu = nan;
    sigma = nan;
    return;
end
fg = fg(:);
bg = bg(:);
M = length(fg);
N = length(bg);
n = M+N;

delta = 1/n;
all = cat(1, bg, fg);
labels = cat(1, bg*0-1, fg*0+1);
[~, od] = sort(all); % ascending
labels = labels(od);
bkpts = find(labels(2:end)-labels(1:end-1));

ai = cat(1, labels(bkpts), -labels(bkpts(end)));
ai(ai>0) = ai(ai>0)*n/M;
ai(ai<0) = ai(ai<0)*(n/N);

% bi is start, ti is end of the i-th section
bi = cat(1, 0, bkpts*delta);
ti = cat(1, bkpts*delta, 1);

Finvbi = norminv(bi);
Finvbi(1) = -1e5;
Finvti = norminv(ti);
Finvti(end) = 1e5;

mu = sum(ai.*(normpdf(Finvbi) - normpdf(Finvti)));
t1=0;
for i=1:length(ai)-1
    aj = ai(i+1:end);
    Finvtj = Finvti(i+1:end);
    Finvbj = Finvbi(i+1:end);
    t1 = t1+ai(i)*...
        sum(aj.*(Finvtj-Finvbj))*...
        (f1(Finvti(i))-f1(Finvbi(i)));
end

t2 = sum(ai.*ai.*Finvti.*(f1(Finvti)-f1(Finvbi)));
t3 = sum(ai.*ai.*(f2(Finvti)-f2(Finvbi)));

A = 2*(t1+t2-t3);
B = (sum(ai.*(f1(Finvti)-f1(Finvbi))))^2;%

sigma = sqrt(A-B)/sqrt(n);

end

function y = f1(x)
    y = x.*normcdf(x)+normpdf(x);
end

function y = f2(x)
    y=0.5*(normcdf(x).*x.^2 - normcdf(x) + normpdf(x).*x);
end
