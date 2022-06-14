function dMap = vid_gradient(vid, smoothness)

if nargin == 1
    smoothness = 1;
end
F = vid * 0;
for i=1:size(F,3)
    F(:,:,i) = imgaussian(vid(:,:,i), smoothness);
end
dx = gradient3(F,'x');
dy = gradient3(F,'y');
dz = gradient3(F,'z');
dMap = sqrt(dx.^2 + dy.^2 + dz.^2);
end