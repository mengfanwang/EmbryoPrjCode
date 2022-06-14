function edge_locs = gapRefine(slabel, idpair, edge_locs)
% given a gap, remove the triangles that outside real gap

[h,w,z] = size(slabel);
edgeMap = false(h,w,z);
edgeMap(edge_locs) = true;

[iy,ix,iz] = ind2sub([h,w,z], edge_locs);
    
z_s = unique(iz);
for i=1:length(z_s)
    curEdgeMap = edgeMap(:,:,z_s(i));
    b = cell(2,1);
    for j=1:2
        tmpStack = slabel(:,:,z_s(i))==idpair(j);
        B = bwboundaries(tmpStack,'noholes');                               % boundary of one labelled region
        % should only one boundary
        if isempty(B)
            error('region removed during the gap testing!');
        end
        yxs = B{1}(1:end-1,:); % the last one is duplicated
        ids = sub2ind([h,w], yxs(:,1), yxs(:,2));                           % index of boundary
        validIds = unique(find(curEdgeMap(ids)>0));                         % overlap of boundary and gap
        if length(validIds)>=2
            yxs = yxs(validIds, :);
            [yy, xx] = find(tmpStack);
            center = [mean(yy) mean(xx)];                                   % center of the labelled region
            dist2c = vecnorm(yxs - center, 2, 2);                           % distance from overlap to region center
            if min(dist2c) > 0 % no need to do any compare
                % scaled to be comparable
                yy2 = (yxs(:,1)-center(1)).*(min(dist2c)./dist2c) + center(1);
                xx2 = (yxs(:,2)-center(2)).*min(dist2c)./dist2c + center(2);
                ym = distanceMat2Mex(yy2, yy2);                             % absolute difference of the two vector
                xm = distanceMat2Mex(xx2, xx2);
                d = sqrt(ym.^2 + xm.^2);

                [b1, b2] = find(d == max(d(:)));

                b{j}(1,:) = yxs(b1(1),:);
                b{j}(2,:) = yxs(b2(1),:);
            end
        end
    end
    if ~isempty(b{1}) && ~isempty(b{2})                                     % change direction 
        if b{1}(1,1) > b{1}(2,1)
            if b{2}(1,1) > b{2}(2,1)
                tmp = b{2}(1,:);
                b{2}(1,:) = b{2}(2,:);
                b{2}(2,:) = tmp;
            end
        end
        b = cat(1, b{:});
        bw = poly2mask(b(:,2), b(:,1),h,w) | slabel(:,:,z_s(i))>0;
        curEdgeMap(~bw) = 0;
    end
    curEdgeMap = imfill(curEdgeMap,'holes');
    edgeMap(:,:,z_s(i)) = curEdgeMap;
end
edge_locs = find(edgeMap);
end