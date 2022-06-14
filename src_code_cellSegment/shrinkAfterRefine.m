function newSeed = shrinkAfterRefine(inLabel, connect, erodeScale)
% using shrinking to find if the fg indeed contains multiple seeds
if nargin == 2
    connect = 6;
    erodeScale = 2;
end
%numSeeds = max(seedRegion(:));
seEff = strel('sphere', erodeScale); % erode on 3d
newLabel = imerode(inLabel, seEff);

% loopCnt = 0;
% while true
%     loopCnt = loopCnt + 1;
%     fg = imerode(fg, seEff);
%     [newSeed, n] = bwlabeln(fg, connect);
%     
%     if loopCnt > 2 || max(seedRegion(newSeed))
%         break;
%     end
%     newSeed(~seedRegion) = 0;
% 
% end

minSzSeed = 10;
stats = regionprops3(newLabel, 'Volume');
if numel()
ismember(newSeed, stats.);
if stats.
newSeed = bwareaopen(newSeed,minSzSeed,connect);


end
