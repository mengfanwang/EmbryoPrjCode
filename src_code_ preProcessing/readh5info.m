function [h5_struct, num_view, name_view] = readh5info(file_name)

    % warning: don't use h5info. It's 200 times slower.
    h5_struct = hdf5info(file_name);
    h5_struct = h5_struct.GroupHierarchy.Groups;
    
    % bdv data format: tTTTTT/sSS/L/cells
    % tTTTTT: time points
    % sSS: id of the setup (view)
    % L: mipmap level (downsample exponent)
    remove_flag = true(size(h5_struct));
    num_view = 0;
    name_view = {};
    for ii = 1:length(h5_struct)
        if h5_struct(ii).Name(2) == 't'
            remove_flag(ii) = false;
        elseif h5_struct(ii).Name(2) == 's'
            num_view = num_view + 1;
            name_view{num_view} = h5_struct(ii).Name(3:end);
        end
    end
    h5_struct(remove_flag) = [];
end