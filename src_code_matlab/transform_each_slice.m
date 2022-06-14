function out3d = transform_each_slice(mat3d)
% transform each slice of a 3d matrix to save the matrix row by row
% 

out3d = zeros(size(mat3d));

for i = 1:size(mat3d,3)
    out3d(:,:, i) = flip(mat3d(:,:,i),1)'; % flip all rows
end

end