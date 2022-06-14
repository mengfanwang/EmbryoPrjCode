function res = transfer_id2cppversion(movieInfo, idx)

ff = isnan(idx);
idx(ff) = 1;
fs = movieInfo.frames(idx);
res = idx;

for i = 1:length(res)
    if fs(i) == 1
        res(i) = res(i) - 1;
    else
        res(i) = res(i) - 1 + sum(movieInfo.n_perframe(1:fs(i)-1));
    end
end


res(ff) = nan;


    