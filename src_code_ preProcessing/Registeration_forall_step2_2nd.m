clc;clear;close all;

data_path = 'H:\Embryo\TM0-49\registration_step1\';
result_path = 'H:\Embryo\TM0-49\registration_step2\';

for tt = 47:49
order = [1 2;2 3; 3 4; 4 1];
tt_ind = num2str(100+tt);
tt_ind = tt_ind(2:3);
tt

diff1 = load([result_path tt_ind '\diff_1.mat']);
diff2 = load([result_path tt_ind '\diff_2.mat']);
diff3 = load([result_path tt_ind '\diff_3.mat']);
diff4 = load([result_path tt_ind '\diff_4.mat']);
y_bias = abs(diff1.y_diff + diff2.y_diff + diff3.y_diff + diff4.y_diff);
x_bias = abs(diff1.x_diff + diff2.x_diff + diff3.x_diff + diff4.x_diff);
z_bias = abs(diff1.z_diff + diff2.z_diff + diff3.z_diff + diff4.z_diff);

% load data
data1 = load([data_path tt_ind '\im_1']);
data2 = load([data_path tt_ind '\im_2']);
data3 = load([data_path tt_ind '\im_3']);
data4 = load([data_path tt_ind '\im_4']);
% add first time registration modification to im1
data2.y_start = data2.y_start + diff1.y_diff;
data2.x_start = data2.x_start + diff1.x_diff;
data2.z_start = data2.z_start + diff1.z_diff;
data2.y_end = data2.y_end + diff1.y_diff;
data2.x_end = data2.x_end + diff1.x_diff;
data2.z_end = data2.z_end + diff1.z_diff;

data3.y_start = data3.y_start + diff1.y_diff + diff2.y_diff;
data3.x_start = data3.x_start + diff1.x_diff + diff2.x_diff;
data3.z_start = data3.z_start + diff1.z_diff + diff2.z_diff;
data3.y_end = data3.y_end + diff1.y_diff + diff2.y_diff;
data3.x_end = data3.x_end + diff1.x_diff + diff2.x_diff;
data3.z_end = data3.z_end + diff1.z_diff + diff2.z_diff;

data4.y_start = data4.y_start + diff1.y_diff + diff2.y_diff + diff3.y_diff;
data4.x_start = data4.x_start + diff1.x_diff + diff2.x_diff + diff3.y_diff;
data4.z_start = data4.z_start + diff1.z_diff + diff2.z_diff + diff3.y_diff;
data4.y_end = data4.y_end + diff1.y_diff + diff2.y_diff + diff3.y_diff;
data4.x_end = data4.x_end + diff1.x_diff + diff2.x_diff + diff3.y_diff;
data4.z_end = data4.z_end + diff1.z_diff + diff2.z_diff + diff3.y_diff;

fprintf('before:');
[diff1.x_diff diff1.y_diff diff1.z_diff;
 diff2.x_diff diff2.y_diff diff2.z_diff;
 diff3.x_diff diff3.y_diff diff3.z_diff;
 diff4.x_diff diff4.y_diff diff4.z_diff;]

pad_size = 10;
data1.im_fuse = padarray(data1.im_fuse, [pad_size pad_size pad_size]);
data2.im_fuse = padarray(data2.im_fuse, [pad_size pad_size pad_size]);
data3.im_fuse = padarray(data3.im_fuse, [pad_size pad_size pad_size]);
data4.im_fuse = padarray(data4.im_fuse, [pad_size pad_size pad_size]);

diff = zeros((2*y_bias+1)^4*(2*x_bias+1)^4*(2*z_bias+1)^4,1);
bias_list = zeros(numel(diff),12); % y im1 im2 im3 im4 / x im1 im2...
bias_count = 0;
for ii12 = -y_bias:y_bias
    for ii23 = -y_bias:y_bias
        for ii34 = -y_bias:y_bias
            for ii41 = -y_bias:y_bias
                for jj12 = -x_bias:x_bias
                    for jj23 = -x_bias:x_bias
                        for jj34 = -x_bias:x_bias
                            for jj41 = -x_bias:x_bias
                                for kk12 = -z_bias:z_bias
                                    for kk23 = -z_bias:z_bias
                                        for kk34 = -z_bias:z_bias
                                            for kk41 = -z_bias:z_bias
if ii12+ii23+ii34+ii41 == -(diff1.y_diff + diff2.y_diff + diff3.y_diff + diff4.y_diff) &&...
   jj12+jj34+jj34+jj41 == -(diff1.x_diff + diff2.x_diff + diff3.x_diff + diff4.x_diff) &&...
   kk12+kk23+kk34+kk41 == -(diff1.z_diff + diff2.z_diff + diff3.z_diff + diff4.z_diff) &&...
   norm([ii12 ii23 ii34 ii41],1) == y_bias && norm([jj12 jj23 jj34 jj41],1) == x_bias && norm([kk12 kk23 kk34 kk41],1) == z_bias
    bias_count = bias_count + 1;
    bias_list(bias_count,:) = [ii12 ii23 ii34 ii41 jj12 jj23 jj34 jj41 kk12 kk23 kk34 kk41];
end
                                            end
                                        end
                                    end
                                end
                            end
                        end
                    end
                end
            end
        end
    end
end
% bias_count = bias_count+1;
fprintf('Count num:%d\n',bias_count);
diff = diff(1:bias_count);
bias_list = bias_list(1:bias_count,:);
if bias_list > 20000
    error('Too much calculation!');
end

tic;
for ind = 1:numel(diff)
    % 12 23 34 14
    bias_temp = bias_list(ind,:);
    diff(ind) = getMSE(data1,data2,pad_size,bias_temp(1:4:9));
    diff(ind) = diff(ind) + getMSE(data2,data3,pad_size,bias_temp(2:4:10));
    diff(ind) = diff(ind) + getMSE(data3,data4,pad_size,bias_temp(3:4:11));
    diff(ind) = diff(ind) + getMSE(data4,data1,pad_size,bias_temp(4:4:12));
end
toc
[~, ind_min] = min(diff);
diff_2nd = bias_list(ind_min,:);
% save([file_path 'diff_2nd'],'diff_2nd');

data1.y_diff = 0;
data1.x_diff = 0;
data1.z_diff = 0;

data2.y_diff = diff1.y_diff + diff_2nd(1);
data2.x_diff = diff1.x_diff + diff_2nd(5);
data2.z_diff = diff1.z_diff + diff_2nd(9);

data3.y_diff = diff1.y_diff + diff2.y_diff + diff_2nd(2);
data3.x_diff = diff1.x_diff + diff2.x_diff + diff_2nd(6);
data3.z_diff = diff1.z_diff + diff2.z_diff + diff_2nd(10);

data4.y_diff = diff1.y_diff + diff2.y_diff + diff3.y_diff + diff_2nd(3);
data4.x_diff = diff1.x_diff + diff2.x_diff + diff3.x_diff + diff_2nd(7);
data4.z_diff = diff1.z_diff + diff2.z_diff + diff3.z_diff + diff_2nd(11);

trans_modification = [data1.x_diff data1.y_diff data1.z_diff;
                      data2.x_diff data2.y_diff data2.z_diff;
                      data3.x_diff data3.y_diff data3.z_diff;
                      data4.x_diff data4.y_diff data4.z_diff;];
fprintf('after');
[diff1.x_diff+diff_2nd(5) diff1.y_diff+diff_2nd(1) diff1.z_diff+diff_2nd(9);
 diff2.x_diff+diff_2nd(6) diff2.y_diff+diff_2nd(2) diff2.z_diff+diff_2nd(10);
 diff3.x_diff+diff_2nd(7) diff3.y_diff+diff_2nd(3) diff3.z_diff+diff_2nd(11);
 diff4.x_diff+diff_2nd(8) diff4.y_diff+diff_2nd(4) diff4.z_diff+diff_2nd(12);]

save([result_path tt_ind '\trans_modification'],'trans_modification');
end

function mse = getMSE(data1,data2,pad_size,diff)
    y_start = max(data1.y_start, data2.y_start + diff(1)) + pad_size;
    y_end   = min(data1.y_end  , data2.y_end   + diff(1)) + pad_size;
    x_start = max(data1.x_start, data2.x_start + diff(2)) + pad_size;
    x_end   = min(data1.x_end  , data2.x_end   + diff(2)) + pad_size;
    z_start = max(data1.z_start, data2.z_start + diff(3)) + pad_size;
    z_end   = min(data1.z_end  , data2.z_end   + diff(3)) + pad_size;
    
    im1 = data1.im_fuse(y_start:y_end,x_start:x_end,z_start:z_end);
    im2 = data2.im_fuse(y_start:y_end,x_start:x_end,z_start:z_end);
    mse = mean((im1-im2).^2,'all');
end