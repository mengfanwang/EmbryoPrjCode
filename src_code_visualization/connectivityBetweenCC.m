function connectTable=connectivityBetweenCC(cellMap)


CC.ImageSize=size(cellMap);
CC.NumObjects=max(cellMap(:));
for cnt=1:CC.NumObjects
    CC.PixelIdxList{cnt}=find(cellMap==cnt);
end

cTable=zeros(CC.NumObjects,3);
szTable=zeros(CC.NumObjects,1);

for cnt=1:CC.NumObjects
    [I1,I2,I3]=ind2sub(CC.ImageSize,CC.PixelIdxList{cnt});
    if length(I1)>1
        cPoint=round(mean([I1,I2,I3]));
    else
        cPoint=round([I1,I2,I3]);
    end
    
    cTable(cnt,:)=cPoint;
    szTable(cnt)=length(CC.PixelIdxList{cnt});
end

D = squareform(pdist(cTable));
sigmaDist=median(szTable)^(1/3);
thres=sigmaDist*4;

connectTable=D<thres;


end