function newIdMap = shapeRefine(idComp, eigCompMap)
% use shape information to decide which gap is real
% INPUT:
% idComp: 3d label map consisting of the components detected by synQuant
% eigCompMap: the principal curvature results (the largest eig value of hessian matrix)
% OUTPUT:
% newIdMap: the label map with all gaps removed

% contact: ccwang@vt.edu; 01/30/2020

[h,w,z] = size(idComp);
newIdMap = zeros(h,w,z);
% figure;
for i=1:z
    slabel = bwlabel(idComp(:,:,i));
    gapMap = eigCompMap(:,:,i);
    olabel = zeros(h,w);
    orgProps = regionprops(slabel, 'Area', 'ConvexArea', 'Orientation', ...
        'MajorAxisLength', 'MinorAxisLength');
    
    for j = 1:numel(orgProps)
        if isempty(orgProps(j).Area)
            continue;
        end
        OrgEllipseScore = 4*orgProps(j).Area / (pi*orgProps(j).MajorAxisLength...
            *orgProps(j).MinorAxisLength);
        if OrgEllipseScore>1
            OrgEllipseScore = 1/OrgEllipseScore;
        end
        OrgCovScore = orgProps(j).Area / orgProps(j).ConvexArea;
        nlabel = slabel==j;%bwlabel(slabel==j);
        nlabel(gapMap) = 0;
        nProps = regionprops(nlabel, 'Area', 'ConvexArea', 'Orientation', ...
            'MajorAxisLength', 'MinorAxisLength');
        if numel(nProps)==0
            olabel(slabel==j) = j;
        else
            es = 0.25*(pi*[nProps.MajorAxisLength].*[nProps.MinorAxisLength]);
            sz = [nProps.Area];
            ElipseScores = min([es;sz])./max([es;sz]);
            covScores = sz./[nProps.ConvexArea];
%             subplot(1,2,1);imshow(slabel==j); title(num2str(min(OrgEllipseScore, OrgCovScore)));
%             subplot(1,2,2);imshow(nlabel); title(num2str(max(ElipseScores)));
            if max(ElipseScores) < OrgEllipseScore || ...
                    (OrgEllipseScore>0.9 && OrgCovScore> 0.85 && orgProps(j).Area < 450)
                olabel(slabel==j) = j;
            else
                olabel(nlabel) = j;
            end
        end
    end
    newIdMap(:,:,i) = olabel;
end
newIdMap = newIdMap>0;
end