function part_ids = partXYpoints(ptsXY, partNum)

ptsX = ptsXY(:,1);
ptsY = ptsXY(:,2);

angles = atan(abs(ptsX)./abs(ptsY));

angles(ptsX>=0 & ptsY<0) = pi - angles(ptsX>=0 & ptsY<0);
angles(ptsX<0 & ptsY<=0) = pi + angles(ptsX<0 & ptsY<=0);
angles(ptsX<=0 & ptsY>0) = 2*pi - angles(ptsX<=0 & ptsY>0);

part_ids = cell(partNum,1);
for i=1:partNum
    part_ids{i} = find(angles>=((i-1)*2*pi/partNum) &...
        angles<(i*2*pi/partNum));
end