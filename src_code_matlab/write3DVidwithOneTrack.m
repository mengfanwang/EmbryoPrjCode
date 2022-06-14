function write3DVidwithOneTrack(orgData,  movieInfo,  trackIdx, p, dataWithParticles, totalClData)
% This function is slightly manipuated for the cell tracking considering
% the lineage for illustartion
%
% write the 3D data with one track and all particles into a given folder
% the particles on the same track are the same color
% INPUT:
% orgData: original 3d video, double but with 0-255
% dataWithParticles: 3d colorful video, the same as orgData
% totalClData: original 3d video data, but with rgb channels

if exist('dataWithParticles','var') == 1
    labelledDataFlag = true;
else
    labelledDataFlag = false;
end
filePath = p.filePath;
fileName = p.fileName;
if ~exist(filePath, 'file')
    mkdir(filePath);
end
[h,w,z,t] = size(orgData);
% remove the ROIs with smaller intensity level
% get the crop of the tracks
leftUp = [h w z];
rightDown = [0 0 0];
fakeFlag = 0;% we remove some overlapped particles (band particles)
cur_track = squeeze(movieInfo.tracks{trackIdx});
if min(size(cur_track))==1 && isempty(find(isnan(cur_track), 1))
    cur_track = cur_track(:);
    parents = [nan;cur_track(1:end-1)];
    cur_track = cat(2, cur_track, parents);
end
for i=1: size(cur_track,1)
    curParticle = cur_track(i,1);
    pt = [movieInfo.yCoord(curParticle), movieInfo.xCoord(curParticle), movieInfo.zCoord(curParticle)];
    if sum(pt<0)>0
        fakeFlag = 1;
    end
    leftUp = min(cat(1,leftUp, pt));
    rightDown = max(cat(1,rightDown, pt));
end
if fakeFlag % particles already been removed because of overlapping
    return;
end

shiftScale = p.shiftScale;%[5 5 3];
leftUp = floor(max(cat(1,leftUp-shiftScale, [1 1 1])));
rightDown = ceil(min(cat(1,rightDown+shiftScale, [h,w,z])));

if exist('totalClData','var')  % we only need to crop the specific region from data in this folder and save
    for i=1:t
        disp(i);
        tmpIm = totalClData{i};
        outIm = tmpIm(leftUp(1):rightDown(1), leftUp(2):rightDown(2), :, leftUp(3):rightDown(3));
        tifwrite(double(outIm)/255, [filePath, '\',fileName, '_', num2str(i)]);
    end
else
    % crop them out
    if ~labelledDataFlag
        orgData = orgData(leftUp(1):rightDown(1), leftUp(2):rightDown(2), leftUp(3):rightDown(3), :);
    else % we have a labelled data with all particles blue
        dataWithParticles = dataWithParticles(leftUp(1):rightDown(1), leftUp(2):rightDown(2), :, ...
            leftUp(3):rightDown(3), :);
    end
    % change the track infor
    for i=1:size(cur_track,1)
        curParticle = cur_track(i,1);
        movieInfo.yCoord(curParticle) = movieInfo.yCoord(curParticle) - leftUp(1) + 1;
        movieInfo.xCoord(curParticle) = movieInfo.xCoord(curParticle) - leftUp(2) + 1;
        movieInfo.zCoord(curParticle) = movieInfo.zCoord(curParticle) - leftUp(3) + 1;
    end
    for i=1:t
        disp(i);
        if ~labelledDataFlag
            orgDataSingle = orgData(:,:,:,i)/p.orgDataDim;
            % build the 3D data
            [h,w, z] = size(orgDataSingle);
            orgData3d = zeros(h,w,3,z);
            for j=1:z
                orgData3d(:,:,1,j) = orgDataSingle(:,:,j);
                orgData3d(:,:,2,j) = orgDataSingle(:,:,j);
                orgData3d(:,:,3,j) = orgDataSingle(:,:,j);
            end
        else
            orgData3d = dataWithParticles(:,:,:,:, i);
            c = size(orgData3d,3);
            if c~=3
                error('dataWithParticles is not colorful data!');
            end
        end
        outIm = imdisplayWithOneTrack3D(orgData3d, movieInfo, trackIdx, i, p);
        tifwrite(outIm/255, fullfile(filePath, [fileName, '_', num2str(i)]));
    end
end