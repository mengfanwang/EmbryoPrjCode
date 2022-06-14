function out_ims = SliceImage(in_im)

[h, w, z] = size(in_im);

h_cut = round( h / 2);
w_cut = round( w / 2 );
if z > 1
    z_cut = round(z/2);
    out_ims = cells(8,1);
    out_ims{1} = in_im(1:h_cut, 1:w_cut, 1:z_cut);
    out_ims{2} = in_im(1:h_cut, 1:w_cut, z_cut+1:end);
    out_ims{3} = in_im(1:h_cut, w_cut+1:end, 1:z_cut);
    out_ims{4} = in_im(1:h_cut, w_cut+1:end, z_cut+1:end);
    out_ims{5} = in_im(h_cut+1:end, 1:w_cut, 1:z_cut);
    out_ims{6} = in_im(h_cut+1:end, 1:w_cut, z_cut+1:end);
    out_ims{7} = in_im(h_cut+1:end, w_cut+1:end, 1:z_cut);
    out_ims{8} = in_im(h_cut+1:end, w_cut+1:end, z_cut+1:end);
else
    out_ims = cells(4,1);
    out_ims{1} = in_im(1:h_cut, 1:w_cut);
    out_ims{2} = in_im(h_cut+1:end, 1:w_cut);
    out_ims{3} = in_im(1:h_cut, w_cut+1:end);
    out_ims{4} = in_im(h_cut+1:end, w_cut+1:end);
end