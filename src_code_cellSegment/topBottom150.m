zz = randn(100000, 500);
zz1 = sort(zz, 2);
% zzb = zz1(:,[1:1500]);
% zzu = zz1(:,[3501:5000]);
cnt = 0;
mu = [];
sigma = [];
a = [];
b = [];
for i=0:0.05:1
    numOv = round(i*size(zz1,2));
    if mod(numOv, 2)>0
        numOv = numOv+1;
    end
    for j = 0.1:0.05:0.9
        div_r = j;
        cnt = cnt+1;
        nonOvNumfg = floor((size(zz1,2)-numOv)*div_r);
        nonOvNumbg = (size(zz1,2)-numOv) - nonOvNumfg;
        rid = nonOvNumbg + randperm(numOv);%sort(nonOvNumbg + randperm(numOv));%
%         rid = sort(rid);%ascending
%         rid = rid(end:-1:1);
        zzb = zz1(:,[1:nonOvNumbg, rid(1:floor(numOv - numOv*div_r))]);
        zzu = zz1(:,[rid(floor(numOv-numOv*div_r)+1:numOv), size(zz1,2)-nonOvNumfg+1:size(zz1,2)]);
        [mu(cnt,1), sigma(cnt,1)] = conventionalIntegralv1(zzb(1,:), zzu(1,:));

%         [order_vec,a(cnt),b(cnt)] = order_info(zzu(1,:), zzb(1,:), rid);
% 
%         [mu(cnt,2), sigma(cnt,2)] = conventionalIntegralv1(zzb(1,:), zzu(1,:), order_vec);
% 
%         [order_vec,a(cnt),b(cnt)] = order_info(zzu(1,:), zzb(1,:));
%         [mu(cnt,3), sigma(cnt,3)] = conventionalIntegralv1(zzb(1,:), zzu(1,:), order_vec);
%         [mu(cnt, 4), sigma(cnt, 4)] = ordStatApprox3sec(zzu(1,:),zzb(1,:),order_vec,a(cnt),b(cnt));
        
        [mu(cnt, 5), sigma(cnt, 5)] = ordStatApproxKsec(zzu(1,:),zzb(1,:));
        zb = mean(zzb,2);
        zu = mean(zzu,2);
        mu(cnt, 4) = mean(zu-zb);
        sigma(cnt, 4) = std(zu-zb);
    end
end
x = 0:0.05:.95;
figure;plot(x, mu(1:end-1,1)); hold on;
plot(x, mu(1:end-1,2)); hold on;
plot(x, mu(1:end-1,3)); hold on;
legend('\mu of L', '\mu of L'' given ov area',...
    '\mu of L'' estimated ov area');
xlabel('overlapping ratio');
ylabel('\mu');

figure;plot(x, sigma(1:end-1,1)); hold on;
plot(x, sigma(1:end-1,2)); hold on;
plot(x, sigma(1:end-1,3)); hold on;
legend('\sigma of L', '\sigma of L'' given ov area',...
    '\sigma of L'' estimated ov area');
xlabel('overlapping ratio');
ylabel('\sigma');

figure;plot(x, b(1:end-1));
hold on; plot(x, x)
xlabel('overlapping ratio');
ylabel('estimated ratio');

x = 0.1:0.05:0.9;
figure;plot(x, mu(1:end,1)); hold on;
plot(x, mu(1:end,2)); hold on;
plot(x, mu(1:end,3)); hold on;
plot(x, mu(1:end,4)); hold on;

legend('\mu of L', '\mu of L'' given ov area',...
    '\mu of L'' estimated ov area');
xlabel('ratio of fg');
ylabel('\mu');

figure;plot(x, sigma(1:end,1)); hold on;
plot(x, sigma(1:end,2)); hold on;
plot(x, sigma(1:end,3)); hold on;
plot(x, sigma(1:end,4)); hold on;
legend('\sigma of L', '\sigma of L'' given ov area',...
    '\sigma of L'' estimated ov area', '3-section approx');
xlabel('ratio of fg');
ylabel('\sigma');

figure;plot(x, b(1:end));
hold on; plot(x, x*0+i)
xlabel('ratio of fg');
ylabel('estimated ratio');


u1 = mean(zzu,2);
b1 = mean(zzb,2);

estVar(1) = var(u1);
estVar(2) = var(b1);
estVar(3) = std(u1-b1);
estVar(4) = var(zzb(randperm(numel(zzb)))'-zzu(:))/size(zzu,2);

zz = randn(10000, 5000);
zz1 = sort(zz, 2);
zzb = zz1(:,[1:150]);
zzu = zz1(:,[351:500]);
rid = 1500+randperm(2000);
zzb = zz1(:, rid(1:1000));
zzu = zz1(:,rid(1001:2000));

u1 = mean(zzu,2);
b1 = mean(zzb,2);

estVar(1) = var(u1);
estVar(2) = var(b1);
estVar(3) = var(u1-b1);

estVar(4) = var(zzb(randperm(numel(zzb)))'-zzu(:))/size(zzu,2);



mu = nan(10000,1);sss = nan(10000,1);

zzu1 = zeros(10000, 1000);
zzb1 = zeros(10000, 1000);
zzu = zeros(10000, 1000);
zzb = zeros(10000, 1000);

for i=1:10000
    % mu(i) = mean(zzu(i,:));
    %[mu(i), sss(i)] = truncatedGauss(0, 1, -inf, zz1(i,150));
    disp(i);
    zz = randn(1, 5000);
    zz1 = sort(zz, 2);
    rid = 1500+randperm(2000);
    zzb(i,:) = zz1(1, rid(1:1000));
    zzu(i,:) = zz1(1,rid(1001:2000));
    lb = zz1(1501);
    ub = zz1(3500);
    t3 = truncate(pd,lb, ub); % lower
    zzu1(i,:) = random(t3, 1,1000);
    zzb1(i,:) = random(t3, 1,1000);
end
zz1 = mean(zzb, 2);
zz2 = mean(zzu, 2);
var(zz1)
var(zz2)
var(zz2-zz1)

zz1 = mean(zzb1, 2);
zz2 = mean(zzu1, 2);


ur = .3;
br = .7;
nn = 1000;
pd = makedist('Normal');
upper_nei = norminv(ur);
lower = norminv(br);
t1 = truncate(pd,-inf, upper_nei); % lower
[n_mu, n_sigma] = truncatedGauss(0, 1, -inf, upper_nei);
t2 = truncate(pd,lower, inf); % upper
[f_mu, f_sigma] = truncatedGauss(0, 1, lower, inf);

t3 = truncate(pd, upper_nei, lower);
z1 = random(t3, 100000,nn);
zz1 = mean(z1, 2);
z2 = random(t3, 100000,nn);
zz2 = mean(z2, 2);

z1 = random(t1, 100000,nn);
zz1 = mean(z1, 2);
z2 = random(t2, 100000,nn);
zz2 = mean(z2, 2);

n_sigma2 = n_sigma*n_sigma/150;
f_sigma2 = f_sigma*f_sigma/150;

realVar(1) = var(zz1);
realVar(2) = var(zz2);
realVar(3) = var(zz2-zz1);
realVar(4) = var(z1(:)-z2(:))/nn;


zzAll = random(pd, 100000, 500);
zzSt = sort(zzAll, 2);
z1 = zzSt(:,1:150);%zzAll;%
%z1(z1>norminv(0.3)) = nan;
zz1 = nanmean(z1, 2);

z1 = zzAll;
z1(z1>norminv(0.3)) = nan;
zz11 = nanmean(z1, 2);
nanNum = sum(~isnan(z1), 2);

z2 = zzSt(:,351:500);%zzAll;%
%z2(z2<norminv(0.7)) = nan;
zz2 = nanmean(z2, 2);

realVar(1) = var(zz1);
realVar(2) = var(zz2);
realVar(3) = var(zz2-zz1);
%realVar(4) = nanvar(z1(~isnan(z1))-z2(~isnan(z2)))/150;

