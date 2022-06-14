function I = display_bnd_map(im, mask)

num_seg = max(mask(:));
b_all = cell(num_seg, 1);
for i=1:num_seg
    b_all{i} = bwboundaries(mask==i);
end
b_all = cat(1, b_all{:});

h = figure; imshow(im); hold on;
for k=1:length(b_all)
   boundary = b_all{k};
   if size(boundary,1) < 10
       continue;
   end
   plot(boundary(:,2), boundary(:,1), 'Color', [0 1 1]/2,'LineWidth',1);
end
hold off;

I = getimage(h);
end