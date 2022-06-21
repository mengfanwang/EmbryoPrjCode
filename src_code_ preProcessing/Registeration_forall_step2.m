clc;clear;close all;

data_path = '\\rs0001\Mengfan\Embryo\TM0-49\registration_step1\';
result_path = '\\rs0001\Mengfan\Embryo\TM0-49\registration_step2\';

for tt = 0:49
order = [1 2;2 3; 3 4; 4 1];
tt_ind = num2str(100+tt);
tt_ind = tt_ind(2:3);
mkdir([result_path tt_ind]);
tt

for oo = 1:4
data1 = load([data_path  tt_ind '\im_' num2str(order(oo,1))]);
data2 = load([data_path  tt_ind '\im_'   num2str(order(oo,2))]);

pad_size = 10;
[x,y,z] = size(data1.im_fuse);
x_start = max(data1.x_start,data2.x_start) + pad_size;
x_end = min(data1.x_end,data2.x_end) + pad_size;
y_start = max(data1.y_start,data2.y_start) + pad_size;
y_end = min(data1.y_end,data2.y_end) + pad_size;
z_start = max(data1.z_start,data2.z_start) + pad_size;
z_end = min(data1.z_end,data2.z_end) + pad_size;

data1.im_fuse = padarray(data1.im_fuse, [pad_size pad_size pad_size]);
data2.im_fuse = padarray(data2.im_fuse, [pad_size pad_size pad_size]);

%
diff = zeros(2*pad_size+1,2*pad_size+1,2*pad_size+1);
tic;
% for x_diff = -pad_size:pad_size
%     fprintf('%d  ',x_diff);
%     for y_diff = -pad_size:pad_size
%         for z_diff = -pad_size:pad_size
%             im1 = data1.im_fuse(y_start:y_end,x_start:x_end,z_start:z_end);
%             im2 = data2.im_fuse(y_start+y_diff:y_end+y_diff,x_start+x_diff:x_end+x_diff,z_start+z_diff:z_end+z_diff);
%             diff(y_diff+pad_size+1,x_diff+pad_size+1,z_diff+pad_size+1) = mean((im1-im2).^2,'all');
%         end
%     end
% end
parfor ind = 1:(2*pad_size+1)^3
    [jj,ii,kk] = ind2sub([21 21 21],ind);
    y_diff = jj - pad_size - 1;
    x_diff = ii - pad_size - 1;
    z_diff = kk - pad_size - 1;
    im1 = data1.im_fuse(y_start:y_end,x_start:x_end,z_start:z_end);
    im2 = data2.im_fuse(y_start+y_diff:y_end+y_diff,x_start+x_diff:x_end+x_diff,z_start+z_diff:z_end+z_diff);
    diff(ind) = mean((im1-im2).^2,'all');
end

[~,min_ind] = min(diff(:));
[y_ind,x_ind,z_ind] = ind2sub(size(diff),min_ind);
y_diff = y_ind - pad_size - 1;
x_diff = x_ind - pad_size - 1;
z_diff = z_ind - pad_size - 1;
fprintf('\n%d %d %d\n', y_diff,x_diff,z_diff);
save([result_path tt_ind '\diff_' num2str(oo)],'diff','y_diff','x_diff','z_diff');

%%
% x_diff = -3; y_diff = 2; z_diff = -3;
color = hsv2rgb([0 1 0.5; 0.25 1 0.5; 0.5 1 0.5; 0.75 1 0.5]);
im1 = zeros(x,y,3,z);
im2 = zeros(x,y,3,z);
for zz = 1:z
    im1(:,:,1,zz) = data1.im_fuse(pad_size+1:end-pad_size,pad_size+1:end-pad_size,pad_size+zz) * color(order(oo,1),1);
    im1(:,:,2,zz) = data1.im_fuse(pad_size+1:end-pad_size,pad_size+1:end-pad_size,pad_size+zz) * color(order(oo,1),2);
    im1(:,:,3,zz) = data1.im_fuse(pad_size+1:end-pad_size,pad_size+1:end-pad_size,pad_size+zz) * color(order(oo,1),3);
    im2(:,:,1,zz) = data2.im_fuse(pad_size+1+y_diff:end-pad_size+y_diff,pad_size+1+x_diff:end-pad_size+x_diff,pad_size+zz+z_diff) * color(order(oo,2),1);
    im2(:,:,2,zz) = data2.im_fuse(pad_size+1+y_diff:end-pad_size+y_diff,pad_size+1+x_diff:end-pad_size+x_diff,pad_size+zz+z_diff) * color(order(oo,2),2);
    im2(:,:,3,zz) = data2.im_fuse(pad_size+1+y_diff:end-pad_size+y_diff,pad_size+1+x_diff:end-pad_size+x_diff,pad_size+zz+z_diff) * color(order(oo,2),3);
end
tifwrite(uint8((im1+im2)/2),[result_path tt_ind '\multiview_xy_' num2str(oo)]);
end
toc
end

% % mean square diff
% x_start = 42; x_end = 294; y_start = 52; y_end = 295;
% win_size = 10;
% diff = zeros(21,21);
% for x_diff = -win_size:win_size
%     for y_diff = -win_size:win_size
%         im1 = data1.im_fuse(x_start:x_end,y_start:y_end,111);
%         im2 = data2.im_fuse(x_start+x_diff:x_end+x_diff,y_start+y_diff:y_end+y_diff,111);
%         diff(x_diff+11,y_diff+11) = mean((im1-im2).^2,'all');
%     end
% end
% 
% 
% im1 = data1.im_fuse(42:294,52:295,111);
% im2 = data2.im_fuse(42-3:294-3,52+2:295+2,111);
% color = hsv2rgb([0 1 0.5; 0.25 1 0.5; 0.5 1 0.5; 0.75 1 0.5]);
% [x,y] = size(im1);
% c1 = zeros(x,y,3);
% c1(:,:,1) = im1*color(2,1);
% c1(:,:,2) = im1*color(2,2);
% c1(:,:,3) = im1*color(2,3);
% c2 = zeros(x,y,3);
% c2(:,:,1) = im2*color(3,1);
% c2(:,:,2) = im2*color(3,2);
% c2(:,:,3) = im2*color(3,3);
% imshow((c1+c2)/300)