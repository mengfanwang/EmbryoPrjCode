function sigma2 = estVarSampleWise(voxVals, OrSt)
% given a set of voxel values, estiamte the vairance they should follow
levels = floor((voxVals-OrSt.stbMinVal)./OrSt.stbUnit) + 1;

if max(levels)>length(OrSt.stbVars)
    error('found super large intensity values');
else
    sigma2 = nanmean(OrSt.stbVars(levels));
end

end