function outIm = imdisplayWithMapColor(orgIm3d, roi3d_org)

V=connectivityBetweenCC(roi3d_org);
IdxLst=mapColor(V);
roi3d=roi3d_org;

for i=1:length(IdxLst)
    roi3d(roi3d_org==i)=IdxLst(i);
end

% generate a 3D colorful data with ROI overlay
[h,w,z] = size(orgIm3d);
% orgIm3d = scale_image(orgIm3d, 0,0.5);
outIm = zeros(h,w,3,z);
colorIntensity = 0.4;

extraColor=round(max(IdxLst)*0.5);
cmap = distinguishable_colors(max(IdxLst)+extraColor);
for cnt=1:extraColor
    [~,i]=max(max(sum(cmap,2),sum(1-cmap,2)));
    cmap(i,:)=[];
end
% cmap=jet(max(IdxLst));
nParticle = max(roi3d(:));

for i=1:z

    tmpROI = roi3d(:,:,i);
    r = zeros(h,w);
    g = zeros(h,w);
    b = zeros(h,w);
    for j=1:nParticle
        r(tmpROI==j) = colorIntensity*cmap(j,1);
        g(tmpROI==j) = colorIntensity*cmap(j,2);
        b(tmpROI==j) = colorIntensity*cmap(j,3);
    end
    outIm(:,:,1,i) = r + orgIm3d(:,:,i);
    outIm(:,:,2,i) = g + orgIm3d(:,:,i);
    outIm(:,:,3,i) = b + orgIm3d(:,:,i);
end
end