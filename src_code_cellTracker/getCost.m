function cost = getCost(movieInfo, cur_id, nei_id)
% get the cost of cur_id --> nei_id
cost = nan;
if movieInfo.frames(cur_id)>movieInfo.frames(nei_id)
    od = find(movieInfo.nei{cur_id} == nei_id,1);
    if lengt(od)~=1
        cost = movieInfo.Cij{cur_id}(od);
    else
        cost = inf;
    end
elseif movieInfo.frames(cur_id)<movieInfo.frames(nei_id)
    od = find(movieInfo.nei{nei_id} == cur_id,1);
    if lengt(od)~=1
        cost = movieInfo.Cij{nei_id}(od);
    else
        cost = inf;
    end
end
end