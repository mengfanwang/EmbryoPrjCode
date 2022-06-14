function [zscore1,zscore2] = regionSig(regMap, vidComp, validNeiMap, OrSt)
%GETORDERSTAT Compute order statistics
% d: transformed data. smask: synapse mask.
% bmask: background mask, MUST be connected. imgVar: image variance

% find balanced groups -----
fg = vidComp(regMap);
[h,w,z] = size(regMap);
% dilate the data (x,y direction) to remove boundary testing
dilate_sc = 2;
vid_di = nan(h+dilate_sc*2, w+dilate_sc*2, z+dilate_sc*2);
vid_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc,...
    dilate_sc+1:end-dilate_sc) = vidComp;
regMap_di = zeros(size(vid_di));
regMap_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc,...
    dilate_sc+1:end-dilate_sc) = regMap;
NeiMap_di = false(size(vid_di));
NeiMap_di(dilate_sc+1:end-dilate_sc, dilate_sc+1:end-dilate_sc,...
    dilate_sc+1:end-dilate_sc) = validNeiMap;
[h,w,z] = size(vid_di);

nei8 = [-h-1, -h, -h+1, -1, 1, h-1, h, h+1];                                                          
%neighbors = [-h*w, nei8, h*w];

nei25 = [-2*h-2:-2*h+2, -h-2, -h+2, -2,2, h-2,h+2, 2*h-2:2*h+2];
% neighbors = [nei8-2*h*w,-2*h*w, nei25-h*w, nei25...
%     , nei25+h*w, 2*h*w, nei8+2*h*w];
neighbors = [nei8-h*w, nei25, nei8+h*w];

surfaceVoxels = find(regMap_di - imerode(regMap_di, strel('cube',3)));
nei_mat = repmat(surfaceVoxels, 1, length(neighbors)) +...
    repmat(neighbors, length(surfaceVoxels), 1);
nei_mat = nei_mat(:);
nei_mat = nei_mat(NeiMap_di(nei_mat));
nei_mat = unique(nei_mat);
fg_neighbors = vid_di(nei_mat);

% approximation of mu and sigma
M = max(length(fg), 10);
N = max(length(fg_neighbors), 10);

[h1, h2] = size(OrSt.NoStbMu);
if (M>=h1 || N>=h2)
    sigmaScl = sqrt((M+N)/500);
    M = floor(M/(M+N)*500);
    N = floor(N/(M+N)*500);
    mu = OrSt.NoStbMu(M, N);
    sigma = OrSt.NoStbSig(M, N);
    sigma = sigma/sigmaScl;
else
    mu = OrSt.NoStbMu(M, N);
    sigma = OrSt.NoStbSig(M, N);
end    

sum_st = mean(fg)-mean(fg_neighbors);
% disp([sum_st mu sigma]);
zscore1 = (sum_st - mu*sqrt(OrSt.NoStbVar))...
    / (sigma*sqrt(OrSt.NoStbVar));

zscore2 = 0;
% [mu, sigma] = conventionalIntegralv1(fg_neighbors, fg);
% disp([mu sigma]);
% zscore2 = (sum_st - mu*sqrt(OrSt.NoStbVar))...
%     / (sigma*sqrt(OrSt.NoStbVar));