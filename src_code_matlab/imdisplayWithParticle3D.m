function outIm = imdisplayWithParticle3D(orgIm3d, movieInfo, curTimePt, p)
% generate a 3D colorful data with particles labelled
[h,w,z] = size(orgIm3d);
particleCl = p.otherParticleCl;
particleSize = p.particleSize;
outIm = zeros(h,w,3,z);
for i=1:z
    ss =orgIm3d(:,:,i);
    outIm(:,:,1,i) = ss;
    outIm(:,:,2,i) = ss;
    outIm(:,:,3,i) = ss;
end
valParticleIdx = find(movieInfo.frames == curTimePt);
for i=1:length(valParticleIdx)
    head = valParticleIdx(i);
    hPt = [movieInfo.yCoord(head),movieInfo.xCoord(head),movieInfo.zCoord(head)];
    outIm = draw3Dparticle(outIm,  hPt, particleSize, particleCl);
end
end