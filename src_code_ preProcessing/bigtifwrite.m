function bigtifwrite(data, file_name)
    % write bigtiff file. for data > 4GB, tifwrite cannot process
    t = Tiff([file_name '.tif'],'w8');
    if size(data,3) == 3
        tagstruct.Photometric = Tiff.Photometric.RGB;
    else
        tagstruct.Photometric = Tiff.Photometric.MinIsBlack;
    end
    tagstruct.Compression = Tiff.Compression.None;
    if isa(data, 'uint8')
        tagstruct.BitsPerSample = 8;
    elseif isa(data, 'uint16')
        tagstruct.BitsPerSample = 16;
    end
    if size(data,4) > 1 && size(data,3) == 3
        tagstruct.SamplesPerPixel = 3;
    else
        tagstruct.SamplesPerPixel = 1;
    end
    tagstruct.SampleFormat = Tiff.SampleFormat.UInt;
    tagstruct.ImageLength = size(data,1);
    tagstruct.ImageWidth = size(data,2);
    tagstruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
    if size(data,4) == 1 % gray image
        for ii = 1:size(data,3)
            t.setTag(tagstruct);
            t.write(uint16(data(:,:,ii)));
            t.writeDirectory();
        end
    elseif size(data,3) == 3
        for ii = 1:size(data,4)
            t.setTag(tagstruct);
            t.write(uint16(data(:,:,:,ii)));
            t.writeDirectory();
        end
    end 
    close(t);