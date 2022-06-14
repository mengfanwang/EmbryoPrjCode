function save_data_in_bin_file(label_maps, threshold_maps, var_maps,... 
eig_maps, save_path, im_file_names)

if nargin == 4
    save_path = ['/home/ccw/Desktop/embryo_res_folder/', ...
        'crop_embryo_data_500x500x30x40/images_downsample'];
    im_file_names = cell(numel(label_maps), 1);
    for i = 1:numel(im_file_names)
        im_file_names{i} = [save_path, '/', num2str(i)];
    end
end


for i = 1:numel(label_maps)
    file_name = [im_file_names{i}, '_%s.bin'];
    % save labels_map as int32 type
    fid = fopen(sprintf(file_name, 'label_map_int32'), 'w'); % open a file, if exist, discard existing content
    fwrite(fid, transform_each_slice(label_maps{i}), 'integer*4');
    fclose(fid);
    
    % save labels_map as uint8 type
    fid = fopen(sprintf(file_name, 'threshold_map_uint8'), 'w'); % open a file, if exist, discard existing content
    fwrite(fid, transform_each_slice(threshold_maps{i}), 'uint8');
    fclose(fid);
    
    % save 2d principal as float type
    fid = fopen(sprintf(file_name, 'principal2d_map_single'), 'w'); % open a file, if exist, discard existing content
    fwrite(fid, transform_each_slice(eig_maps{i}{1}), 'float32');
    fclose(fid);
    
    % save 3d principal as float type
    fid = fopen(sprintf(file_name, 'principal3d_map_single'), 'w'); % open a file, if exist, discard existing content
    fwrite(fid, transform_each_slice(eig_maps{i}{2}), 'float32');
    fclose(fid);
    
    % save variance map as float type
    fid = fopen(sprintf(file_name, 'var_map_single'), 'w'); % open a file, if exist, discard existing content
    fwrite(fid, transform_each_slice(var_maps{i}{1,1}), 'float32');
    fclose(fid);
    
    % save stablized variance map as float type
    fid = fopen(sprintf(file_name, 'stb_var_map_single'), 'w'); % open a file, if exist, discard existing content
    fwrite(fid, transform_each_slice(var_maps{i}{1,2}), 'float32');
    fclose(fid);
    
    % save var and var trend as float type
    fid = fopen(sprintf(file_name, 'var_trend_single'), 'w'); % open a file, if exist, discard existing content
    fwrite(fid, var_maps{i}{2,1}, 'float32');
    fwrite(fid, var_maps{i}{3,1}, 'float32');
    fclose(fid);
    
    % save save stablized var and var trend as float type as float type
    fid = fopen(sprintf(file_name, 'stb_var_trend_single'), 'w'); % open a file, if exist, discard existing content
    fwrite(fid, var_maps{i}{2,2}, 'float32');
    fwrite(fid, var_maps{i}{3,2}, 'float32');
    fclose(fid);
end
end