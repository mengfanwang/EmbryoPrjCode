function outData = generate3DImWithParticle(orgData, movieInfo, p)
% write the 3D data into a given folder
[h,w,z,t] = size(orgData);
outData = zeros(h,w, 3,z,t);
% remove the ROIs with smaller intensity level
orgDataDim = p.orgDataDim;
for i=1:t
    disp(i);
    orgDataSingle = orgData(:,:,:,i)/orgDataDim;
    outIm = imdisplayWithParticle3D(orgDataSingle, movieInfo, i, p);
    outData(:,:,:,:, i) = outIm;
end
end