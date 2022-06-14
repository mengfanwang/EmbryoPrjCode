clc;clear;close all;
dbstop if error
addpath('..\src_code_cellSegment');


%% simulation
mean_all = zeros(1,50);
std_all = zeros(1,50);
org_std_all = zeros(1,50);
for ii = 1:50
    mm = ii/20;
len = 100;
% assume image follow gaussian distribution 
    vid = normrnd(0, 1, [len+4,len+4, len+4]);
    sigma = 0;
    [eig_all, ~] = principalCv3d(vid, zeros(size(vid)), mm, ones(size(vid)));
    eig_all = eig_all(3:end-2, 3:end-2, 3:end-2);
% assume principal curvature follow gaussian
%     [eig_all, org_all] = priCvt3d_from_randGrad(len, len, len, mm, 'norm');
    
    
mean_all(ii) = mean(eig_all(:));
std_all(ii) = std(eig_all(:));
% org_std_all(ii) = std(org_all(:));
% histogram(eig_all(:));hold on;

% a = normrnd(mean(eig_all(:)), std(eig_all(:)), len^3, 1);
% histogram(a(:));
end
plot(org_std_all, mean_all); hold on;
plot(org_std_all, std_all); hold on;    


%% symbol calculation
% gradient calculation method: Gi = (X_(i+1) - X_(i-1))/2
X = sym('X%d%d', [5 5 5]);
Gx = gradient_x(X);
Gy = gradient_y(X);
Gz = gradient_z(X);
Gxx = gradient_x(Gx);
Gxy = gradient_y(Gx);
Gxz = gradient_z(Gx);
Gyy = gradient_y(Gy);
Gyz = gradient_z(Gy);
Gzz = gradient_z(Gz);
H = [Gxx Gxy Gxz;...
     Gxy Gyy Gyz;...
     Gxz Gyz Gzz;];
 
syms lambda
f = det(H - lambda*eye(3,3));
eig = solve(f==0,lambda);

%% eigenvalue calculation
syms a b c d e f lambda
syms A B C D
H = [a b c; b d e; c e f;] - eye(3,3)*lambda;
coef = coeffs(det(H),lambda);
A = coef(4);
B = coef(3);
C = coef(2);
D = coef(1);
p1 = B*C/6 + B^3/27 + D/2;
p2 = C/3 + B^2/9;

%% pdf cakculation


function Gx = gradient_x(X)
    Gx = (X(3:end, 2:end-1, 2:end-1) - X(1:end-2, 2:end-1, 2:end-1))/2;
end

function Gy = gradient_y(X)
    Gy = (X(2:end-1, 3:end, 2:end-1) - X(2:end-1, 1:end-2, 2:end-1))/2;
end

function Gz = gradient_z(X)
    Gz = (X(2:end-1, 2:end-1, 3:end) - X(2:end-1, 2:end-1, 1:end-2))/2;
end

function [eig_all, org_all] = priCvt3d_from_randGrad(h,w,z,sigma, distribution)

if strcmp(distribution, 'norm')
    Dxx = normrnd(0, sigma, [h, w, z]);
    Dxy = normrnd(0, sigma, [h, w, z]);
    Dyy = normrnd(0, sigma, [h, w, z]);
    Dxz = normrnd(0, sigma, [h, w, z]);
    Dyz = normrnd(0, sigma, [h, w, z]);
    Dzz = normrnd(0, sigma, [h, w, z]);
elseif strcmp(distribution, 'poiss')
    Dxx = poissrnd(sigma, [h, w, z]) - sigma;
    Dxy = poissrnd(sigma, [h, w, z]) - sigma;
    Dyy = poissrnd(sigma, [h, w, z]) - sigma;
    Dxz = poissrnd(sigma, [h, w, z]) - sigma;
    Dyz = poissrnd(sigma, [h, w, z]) - sigma;
    Dzz = poissrnd(sigma, [h, w, z]) - sigma;
elseif strcmp(distribution, 'uniform')
    Dxx = (rand([h, w, z]) - 0.5) * sigma;
    Dxy = (rand([h, w, z]) - 0.5) * sigma;
    Dyy = (rand([h, w, z]) - 0.5) * sigma;
    Dxz = (rand([h, w, z]) - 0.5) * sigma;
    Dyz = (rand([h, w, z]) - 0.5) * sigma;
    Dzz = (rand([h, w, z]) - 0.5) * sigma;
elseif strcmp(distribution, 'trunNorm')
    Dxx = normrnd(0, sigma, [h, w, z]);
    Dxy = normrnd(0, sigma, [h, w, z]);
    Dyy = normrnd(0, sigma, [h, w, z]);
    Dxz = normrnd(0, sigma, [h, w, z]);
    Dyz = normrnd(0, sigma, [h, w, z]);
    Dzz = normrnd(0, sigma, [h, w, z]);
    Dxx(Dxx<0) = 0;
    Dxy(Dxy<0) = 0;
    Dyy(Dyy<0) = 0;
    Dyz(Dyz<0) = 0;
    Dxz(Dxz<0) = 0;
    Dzz(Dzz<0) = 0;
elseif strcmp(distribution, 'GOE')
    Dxx = sqrt(2) * normrnd(0, sigma, [h, w, z]);
    Dxy = normrnd(0, sigma, [h, w, z]);
    Dyy = sqrt(2) * normrnd(0, sigma, [h, w, z]);
    Dxz = normrnd(0, sigma, [h, w, z]);
    Dyz = normrnd(0, sigma, [h, w, z]);
    Dzz = sqrt(2) * normrnd(0, sigma, [h, w, z]);
end


% test each connected component
eig_all = zeros(h,w,z);
org_all = zeros(h,w,z);
eig1 = zeros(h,w,z);

dir_y = zeros(h,w,z);
dir_x = zeros(h,w,z);
dir_z = zeros(h,w,z);

fmap = ones(h,w,z);

s = regionprops3(fmap, {'VoxelIdxList'});
for i=1:numel(s.VoxelIdxList)%[3 47 76]%
    %disp(i);
    vox = s.VoxelIdxList{i};
    xx = Dxx(vox); yy = Dyy(vox); zz = Dzz(vox);
    xy = Dxy(vox); xz = Dxz(vox); yz = Dyz(vox);
    
    C = zeros(numel(s.VoxelIdxList{i}),3);
    dir_xyz = zeros(numel(s.VoxelIdxList{i}),3);
    parfor j=1:numel(s.VoxelIdxList{i})
        MM = [xx(j), xy(j), xz(j);...
            xy(j), yy(j), yz(j);...
            xz(j), yz(j), zz(j)];
        [Evec,Eval] = eig(MM);
        dEval = diag(Eval);
        %[~,od] = sort(abs(dEval),'descend');
        %C(j,:) = dEval(od)';
        [c,od] = sort(dEval,'descend');
        C(j,:) = c';
        dir_xyz(j,:) = Evec(:, od(1))';
    end
    dir_x(vox) = dir_xyz(:,1);
    dir_y(vox) = dir_xyz(:,2);
    dir_z(vox) = dir_xyz(:,3);
    eig_all(vox) = C(:,1);
    org_all(vox) = Dxx;
end
end