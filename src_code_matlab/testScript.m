
zz = randn(100000, 500);

zz1 = sort(zz, 2);

zzb = zz1(:,[1:200, 251:300]);
zzu = zz1(:,[201:250, 301:500]);

rid = 150+randperm(200);
zzb = zz1(:,[1:150, rid(1:100)]);
zzu = zz1(:,[rid(101:200), 351:500]);

zzb = zz1(:,[1:150]);
zzu = zz1(:,[351:500]);

kk = mean(zzu,2)-mean(zzb,2);
var(zzu(:))/150
var(zzb(:))/150
var(mean(zzb,2))
var(mean(zzu,2))
var(kk)

var(zzb(randperm(numel(zzb)))'-zzu(:))/size(zzu,2)


zzb = zz1(:,[1:150, 251:350]);
zzu = zz1(:,[151:250, 351:500]);

zzb = zz1(:,1:150);
zzu = zz1(:,351:500);

zzb = zz1(:,[1:80, 101:120]);
zzu = zz1(:,[81:100, 121:500]);

zzb = zz1(:,1:2:500);
zzu = zz1(:,2:2:500);

%% 
kk = randn(100000, 1);
kk1 = kk(kk<norminv(0.3));
kk2 = kk(kk<norminv(0.7) & kk>norminv(0.3));
kk3 = kk2(1:round(length(kk2)/2));
kk0 = [kk1; kk3];
sb = var(kk0)/250;

kk4 = kk(kk>norminv(0.7));
kk5 = kk2(round(length(kk2)/2):end);
kk6 = [kk4; kk5];
su = var(kk6)/250;

ss0 = sqrt((sb+su)*0.6);

ss = std(mean(zzu, 2)-mean(zzb, 2))

[~, sigma] = conventionalIntegralv1(zzb(1,:), zzu(1,:))


%% ovTruncated Gaussian
N = size(zzb,2);
M = size(zzu,2);
fg = zzu(1,:)';
fg_neighbors = zzb(1,:)';
all_v = cat(1, fg, fg_neighbors);
labels = cat(1, ones(M, 1), 2*ones(N, 1));
[~, od_v] = sort(all_v, 'descend');
labels = labels(od_v);
ratio_v = (1:length(labels))'./length(labels);
fg_ratio = mean(ratio_v(labels==1));
lower_fg = -norminv(2*fg_ratio);
[f_mu, f_sigma] = truncatedGauss(0, 1, lower_fg, inf);
f_sigma = f_sigma/sqrt(M);

nei_ratio = 2*(1-mean(ratio_v(labels==2)));
upper_nei = norminv(nei_ratio);
[n_mu, n_sigma] = truncatedGauss(0, 1, -inf, upper_nei);
n_sigma = n_sigma/sqrt(N);

mu = f_mu - n_mu;
sigma = sqrt(f_sigma^2 + n_sigma^2)

%%
[~, s1] = truncatedGauss(0, 1, -inf, 0.3);
[~, s2] = truncatedGauss(0, 1, 0.7, inf);
a = sqrt(s1^2/N+s2^2/M);

[~, s1] = truncatedGauss(0, 1, -inf, norminv(0.30));
[~, s2] = truncatedGauss(0, 1, norminv(1-0.40), inf);
a = sqrt(s1^2*0.8/N+s2^2*0.8/M);
[~, s3] = truncatedGauss(0, 1, norminv(0.40), norminv(1-0.40));
b = sqrt(s3^2*0.2/N+s3^2*0.2/M);
s4 = sqrt(a^2+b^2)

[~, s1] = truncatedGauss(0, 1, -inf, norminv(0.52));
[m1, s1] = truncatedGauss(0, 1, -inf, norminv(0.3))
[~, s1] = truncatedGauss(0, 1, norminv(0.40), norminv(0.60))
[~, s2] = truncatedGauss(0, 1, norminv(1-0.40), inf);
a = sqrt(s1^2*0.8/N+s2^2*0.8/M);




pd = makedist('Normal');
upper_nei = norminv(0.3);
lower = norminv(0.7);
t1 = truncate(pd,-inf, upper_nei); % lower
[n_mu, n_sigma] = truncatedGauss(0, 1, -inf, upper_nei);
t2 = truncate(pd,lower, inf); % upper
[f_mu, f_sigma] = truncatedGauss(0, 1, lower, inf);

z1 = random(t1, 100000,150);
zz1 = mean(z1, 2);
z2 = random(t2, 100000,150);
zz2 = mean(z2, 2);

n_sigma2 = n_sigma*n_sigma/150;
f_sigma2 = f_sigma*f_sigma/150;

realVar(1) = var(zz1);
realVar(2) = var(zz2);
realVar(3) = var(zz2-zz1);
realVar(4) = var(z1(:)-z2(:))/150;
[~, sigma] = conventionalIntegralv1(z1(1,:), z2(1,:))
%% 

t11 = truncate(pd,-inf, norminv(0.4)); % lower
t12 = truncate(pd,0, norminv(0.6)); % lower

t21 = truncate(pd,norminv(0.6), inf); % upper
t22 = truncate(pd,norminv(0.4), 0); % upper


z11 = random(t11, 100000,200);
z12 = random(t12, 100000,50);
z1 = [z11, z12];
zz1 = mean(z1, 2);
z21 = random(t21, 100000,200);
z22 = random(t22, 100000,50);
z2 = [z21, z22];
zz2 = mean(z2, 2);

realVar(1) = var(zz1);
realVar(2) = var(zz2);
realVar(3) = var(zz2-zz1);

[~, sigma] = conventionalIntegralv1(z1(1,:), z2(1,:));
realVar(4) = sigma^2;
