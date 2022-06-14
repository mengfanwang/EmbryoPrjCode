function [sigma2, xx] = temporalVar(vid, vid2, varRatio)

vid_stb = sqrt(vid+3/8);

se = ones(5,5,5);
se(2:4,2:4,2:4) = 0;
se(:,:,1:2) = 0;
se(:,:,4:5) = 0;
se = se./sum(se(:));

diff = sqrt(vid2+3/8)-sqrt(vid+3/8);
vid_stb_s = imboxfilt3(vid_stb,[3,3,3]);%imfilter(vid_stb,se);%vid_stb;%

% figure;histogram(diff(vid_stb_s>a & vid_stb_s<b))
% var_t = var(diff(vid_stb_s>a & vid_stb_s<b))/2;

levels = 200;
unit = (max(vid_stb_s(:)) - min(vid_stb_s(:)))/(levels-1);
varDifferentInt = nan;
localmean2x = floor(vid_stb_s ./ unit) + 1;

eleCnt = 0;
ele0 = length(find(localmean2x == 1));
xx = zeros(levels,1);
for i = 2: levels
    %fprintf([num2str(i) '/' num2str(levels) '\n'])
    cur_loc = find(localmean2x == i);
    varDifferentInt = diff(cur_loc);
    xx(i) = var(varDifferentInt)/2;
    eleCnt = eleCnt + length(cur_loc);
    if eleCnt / (numel(localmean2x) - ele0) > varRatio
        break;
    end
end
sigma2 = var(varDifferentInt);


end