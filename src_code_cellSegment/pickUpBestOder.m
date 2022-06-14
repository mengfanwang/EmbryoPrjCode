function od = pickUpBestOder(z, z_debias, sz_nei)
% given zscore, debiased tscore, and neighboring size, pick up the best
% threshold 

[~, max_loc] = nanmax(sz_nei);

[max_debias, max_debias_loc] = nanmax(z_debias(1:max_loc));

[~, max_z_loc] = nanmax(z(max_loc:end));

max_z_loc = max_z_loc + max_loc - 1;

if z_debias(max_z_loc) > max_debias
    od = max_z_loc;
else
    od = max_debias_loc;
end
end