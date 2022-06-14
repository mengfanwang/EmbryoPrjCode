function orgIm3d = imdisplayWithOneTrack3D(orgIm3d, movieInfo, trackIdx, timePt, p)
% generate a 3D colorful data with only one Track
% can be merged to imdisplayWithTrack3D
maxInt = p.maxInt;%max(orgIm3d(:));
particleSize = p.particleSize;
particleCl = p.particleCl ;
stoppedParticleCl = p.stoppedParticleCl;
lineWidth = p.lineWidth;
RegLineCl = p.RegLineCl;
stoppeLineCl = p.stoppeLineCl;

cur_track = squeeze(movieInfo.tracks{trackIdx});
if min(size(cur_track))==1 && isempty(find(isnan(cur_track), 1))
    cur_track = cur_track(:);
    parents = [nan;cur_track(1:end-1)];
    cur_track = cat(2, cur_track, parents);
end

for i=trackIdx % this is the only difference from function imdisplayWithTrack3D
    fms = movieInfo.frames(cur_track(:,1));
    % continuous path
    node_loc = find(fms==timePt);% if > 1 means seperation
    if ~isempty(node_loc) 
        if size(cur_track,1)==1 % the track consist of only one particle
            head = cur_track(1,1);
            pt = [movieInfo.yCoord(head), movieInfo.xCoord(head), movieInfo.zCoord(head)];
            orgIm3d = draw3Dparticle(orgIm3d,  pt, particleSize, particleCl);
        elseif timePt==1 || node_loc(1)==1 % start of the track
            for j=1:length(node_loc)
                head = cur_track(node_loc(j),1);
                hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
                orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, particleCl);
            end
        else %
            node_loc = node_loc(end);
            for j = node_loc:-1:2
                tail = cur_track(j,2);
                head = cur_track(j,1);
                tPt = [movieInfo.yCoord(tail),movieInfo.xCoord(tail), movieInfo.zCoord(tail)];
                hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
                orgIm3d = draw3Dparticle(orgIm3d,  tPt, particleSize, particleCl);
                orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, particleCl);
                orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, RegLineCl);
            end
        end
    else
        % if the path is jumping this frame
        if sum(fms>timePt) && sum(fms<timePt) % make sure not stop
            node_loc = find(fms<timePt);
            node_loc = node_loc(end);
            if node_loc==1
                head = cur_track(1,1);
                hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
                orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
            else
                for j = node_loc:-1:2
                    tail = cur_track(j,2);
                    head = cur_track(j,1);
                    tPt = [movieInfo.yCoord(tail),movieInfo.xCoord(tail), movieInfo.zCoord(tail)];
                    hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
                    orgIm3d = draw3Dparticle(orgIm3d,  tPt, particleSize, stoppedParticleCl);
                    orgIm3d = draw3Dparticle(orgIm3d,  hPt, particleSize, stoppedParticleCl);
                    orgIm3d = draw3Dline(orgIm3d, hPt, tPt, lineWidth, stoppeLineCl);
                end
            end
        end
    end
end

end