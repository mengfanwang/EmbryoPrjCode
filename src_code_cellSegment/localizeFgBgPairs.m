function [pairs, pairs_locs] = localizeFgBgPairs(vid, regMap, validNei, ...
    bandWidth, minSz, pairNum)
% fgMap/bgMap: indicating foreground/background pixels
% minSz: minimum size in each partition
% pairNum: number of partitions we want. E.g. 4, 8

if nargin<6
    pairNum=8;
end

[h,w,slice] = size(vid);

sph = strel('sphere', bandWidth);
se = strel(sph.Neighborhood(:,:,bandWidth+1));
nei_locs = imdilate(regMap, se)>regMap & validNei; % h*w*slice
fg_locs = regMap & imdilate(nei_locs, se);

pairs = cell(100*pairNum*slice,2);
pairs_locs = cell(100*pairNum*slice,2);
pair_cnt = 0;
for i=1:slice
    cur_slice_vid = vid(:,:,i);
    cur_slice_fg = fg_locs(:,:,i);
    cur_slice_nei = nei_locs(:,:,i);
    [l, n] = bwlabel(cur_slice_fg,8);
    s = regionprops(l,'PixelList');
    zz = zeros(h,w);
    zz(regMap(:,:,i)>cur_slice_fg)=1;
    for j=1:n
        if n>1
            fg_pixXY = s(j).PixelList;
        else
            fg_pixXY = s.PixelList;
        end
        centerXY = mean(fg_pixXY);
        fg_parts = partXYpoints(fg_pixXY-centerXY, pairNum);
        fg_len = cellfun(@length, fg_parts);
        [yy,xx] =  find(imdilate(l==j, se) & cur_slice_nei);
        nei_pixXY = [xx, yy];
        nei_parts = partXYpoints(nei_pixXY-centerXY, pairNum);
        nei_len = cellfun(@length, nei_parts);
        
        validPair = find(min([fg_len nei_len], [], 2) >= minSz);
        
        if isempty(validPair) && min(sum(fg_len), sum(nei_len))>minSz
            pair_cnt = pair_cnt + 1;
            tmp_fg_locs = sub2ind([h,w], ...
                fg_pixXY(:,2), fg_pixXY(:,1));
            pairs{pair_cnt,1} = cur_slice_vid(tmp_fg_locs);
            pairs_locs{pair_cnt,1} = tmp_fg_locs + (i-1)*h*w;
            
            tmp_bg_locs = sub2ind([h,w], ...
                nei_pixXY(:,2), nei_pixXY(:,1));
            pairs{pair_cnt,2} = cur_slice_vid(tmp_bg_locs);
            pairs_locs{pair_cnt,2} = tmp_bg_locs + (i-1)*h*w;
        else
            for k = 1:length(validPair)
                pair_cnt = pair_cnt + 1;
                tmp_fg_locs = sub2ind([h,w], ...
                    fg_pixXY(fg_parts{validPair(k)},2), ...
                    fg_pixXY(fg_parts{validPair(k)},1));
                pairs{pair_cnt,1} = cur_slice_vid(tmp_fg_locs);
                pairs_locs{pair_cnt,1} = tmp_fg_locs + (i-1)*h*w;
                
                tmp_bg_locs = sub2ind([h,w], ...
                    nei_pixXY(nei_parts{validPair(k)},2), ...
                    nei_pixXY(nei_parts{validPair(k)},1));
                pairs{pair_cnt,2} = cur_slice_vid(tmp_bg_locs);
                pairs_locs{pair_cnt,2} = tmp_bg_locs + (i-1)*h*w;
                zz(tmp_fg_locs) = k*2;
                zz(tmp_bg_locs) = k*2+1;
            end
            figure;imshow(label2rgb(zz));
        end
    end
end

pairs = pairs(1:pair_cnt,:);
pairs_locs = pairs_locs(1:pair_cnt,:);
end