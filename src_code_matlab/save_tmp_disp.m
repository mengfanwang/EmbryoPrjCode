function save_tmp_disp(movieInfo, org_refine_res, tr_id, dir_folder)
% save cropped region consisting of segmentation results
%INPUT:
% movieInfo:
% org_refine_res:
% tr_id: track id
% nd_id: node id
% time_sft: show frames in [t-time_sft:t+time_sft]
%OUTPUT:

if nargin < 4
    dir_folder = [];
end
[h,w,z] = size(org_refine_res{1});
leftUp = [h w z];
rightDown = [0 0 0];
fakeFlag = 0;% we remove some overlapped particles (band particles)
for i=1: length(movieInfo.tracks{tr_id})
    curParticle = movieInfo.tracks{tr_id}(i);
    if isnan(curParticle) || curParticle < 1
        continue;
    end
    pt = [movieInfo.yCoord(curParticle), movieInfo.xCoord(curParticle), movieInfo.zCoord(curParticle)];
    if sum(pt<0)>0
        fakeFlag = 1;
    end
    leftUp = min(cat(1,leftUp, pt));
    rightDown = max(cat(1,rightDown, pt));
end
shiftScale = [25 25 10];
leftUp = floor(max(cat(1,leftUp-shiftScale, [1 1 1])));
rightDown = ceil(min(cat(1,rightDown+shiftScale, [h,w,z])));
cur_frs = movieInfo.frames(movieInfo.tracks{tr_id});
if exist(fullfile(dir_folder, 'tmp_disp'),'dir')
    mkdir(fullfile(dir_folder, 'tmp_disp'));
end
for i = min(cur_frs):max(cur_frs)
    cut_im = org_refine_res{i}(leftUp(1):rightDown(1), leftUp(2):rightDown(2), leftUp(3):rightDown(3));
    fr1 = label2rgb3d(cut_im, 'jet', [0 0 0], 'shuffle');
    tifwrite(fr1, fullfile(dir_folder, 'tmp_disp', [num2str(tr_id),'_', num2str(i)]));
end

end