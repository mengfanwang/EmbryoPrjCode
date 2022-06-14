% k1 embro data
%output_file = 'output_embryo_k1.txt';
% k2 embro data
output_file = 'output_embryo_k2.txt';
track_vec = dlmread(output_file);

addpath('/home/ccw/Dropbox/cc_ImHandle/');

locs = find(track_vec<1);
if length(track_vec) > locs(end)
    track_vec = track_vec(1:locs(end));
end

track_len = locs(2:end) - locs(1:end-1) - 1;
track_len = [locs(1)-1; track_len]/2;
track_vec(locs) = [];
track_vec = track_vec(1:2:end)/2;
tr_start_pts = track_vec([0; track_len(1:end-1)] + 1);

trajectories = mat2cell(track_vec, track_len,1);



load('embryo_10tb_detections.mat');
for i=1:numel(detections)
    detections{i} = cat(2, detections{i}, i*ones(size(detections{i},1),1));
end
det = cat(1, detections{5:end});

det_num_per_frame = cellfun(@length, detections);

numTracks = numel(trajectories);
l_tracks = cellfun(@length, trajectories);

st_1st = [];%find(tr_start_pts < det_num_per_frame(5));
for i=1:numTracks
    if trajectories{i}(1) < det_num_per_frame(5)
        st_1st = [st_1st; i];
    end
end


[~, od] = sort(l_tracks(st_1st),'descend');
track_display = trajectories(st_1st(od(1:200)));

% time_occupy = nan(numel(track_display), 600);
% for i=1:numel(track_display)
%     tmp = det(track_display{i},4)';
%     time_occupy(i,1:length(tmp)) = tmp;
% end



%det(:,3) = ceil(det(:,3)/100);


%% read image
det_org = cat(1, detections{5:end});
%det(:,1:3) = ceil(det(:,1:3)./5);
det_org(:,1) = det_org(:,1) - 759;
%det(:,2) = det(:,2) - 759;
det_org(:,3) = det_org(:,3) - 268;
luts;
inValid = [];
for i=1:numel(track_display)
    if ~isempty(find(det(track_display{i},1)>1301 | det(track_display{i},3)>419, 1))
        inValid = [inValid, i];
    end
end
track_display(inValid) = [];
track_display = track_display(1:100);

det = det_org;
det(:,1:2) = ceil(det(:,1:2)./5);
p = [];
p.particleSize = 2;
p.lineWidth = 1;
track_display = track_display(1);
p.cmap = hsv(numel(track_display));
p.cmap = uint8(p.cmap*255);

tt = det(track_display{1},end);
for ff = 1:50:length(tt)
    f = tt(ff);
    disp(f);
    im_path = sprintf('/media/congchao/New Volume/ImageData/Mmu_E1_CAGTAG1.TM000%03d_timeFused_blending/',f);
    %im_name = 'SPM00_TM000082_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb';
    im_name = sprintf('SPM00_TM000%03d_CM00_CM01_CHN00.fusedStack.corrected.shifted.klb',f);
    
    [im, header] = readKLBstack(fullfile(im_path, im_name), 4);
    
    im = double(im(760:2060,:,269:687));
    im = double(im)/double(max(im(:)));
    tifwrite(im,sprintf('track_dispaly_%03d_org_big',f));
    continue;
    [~,~,z] = size(im);
    [h,w] = size(imresize(im(:,:,1),0.2));
    im_rz = zeros(h,w,z);
    for i=1:z
        im_rz(:,:,i) = imresize(im(:,:,i),0.2);
    end
    im_fg = uint8(255*im_rz/max(im_rz(:)));
    %tifwrite(im_fg,sprintf('track_dispaly_%03d_org',f));
    [h,w,z] = size(im_fg);
    orgIm3d = uint8(zeros(h,w,3,z));
    for i=1:z
        disp(i);
        cur_im = im_fg(:,:,i);

        tmp = uint8(zeros(h,w));
        tmp(:) = royal(cur_im(:)+1,2);
        orgIm3d(:,:,1,i) = tmp;

        tmp(:) = royal(cur_im(:)+1,3);
        orgIm3d(:,:,2,i) = tmp;

        tmp = uint8(zeros(h,w));
        tmp(:) = royal(cur_im(:)+1,4);
        orgIm3d(:,:,3,i) = tmp;
    end
    %tifwrite(orgIm3d,sprintf('track_dispaly_%03d_blue',f));

%     im_rz = uint8(1.5*255*im_rz/max(im_rz(:)));
%     [h,w,z] = size(im_rz);
%     orgIm3d = uint8(zeros([h,w,3,z]));
%     for i=1:size(orgIm3d,4)
%         orgIm3d(:,:,1,i) = im_rz(:,:,i);
%         orgIm3d(:,:,2,i) = im_rz(:,:,i);
%         orgIm3d(:,:,3,i) = im_rz(:,:,i);
%     end
    
%     outIm = drawTrack_klb(orgIm3d, track_display, det, p, f);
    pt = det(track_display{1}(tt==f),1:3);
    orgIm3d = draw3Dparticle(orgIm3d,  pt, p.particleSize, p.cmap(1,:));
    tifwrite(orgIm3d,sprintf('track_dispaly_%03d',f));

end