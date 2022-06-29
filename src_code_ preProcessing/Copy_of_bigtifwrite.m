function bigtifwrite(data, file_name)
    % write bigtiff file. for data > 4GB, tifwrite cannot process
    t = Tiff(file_name,'w8');
    if size(data,3) == 3
        setTag(t,'Photometric',Tiff.Photometric.RGB);
    else
        setTag(t,'Photometric',Tiff.Photometric.MinIsBlack);
    end
    setTag(t,'Compression',Tiff.Compression.None);
    if isa(data, 'uint8')
        setTag(t,'BitsPerSample',8);
    elseif isa(data, 'uint16')
        setTag(t,'BitsPerSample',16);
    end
%     if size(data,3) == 3
%         setTag(t,'SamplesPerPixel',3);
%     else
%         setTag(t,'SamplesPerPixel',1);
%     end
    setTag(t,'SamplesPerPixel',size(data,3));
    setTag(t,'SampleFormat',Tiff.SampleFormat.UInt);
    setTag(t,'ExtraSamples',Tiff.ExtraSamples.Unspecified);
    setTag(t,'ImageLength',size(data,1));
    setTag(t,'ImageWidth',size(data,2));
    setTag(t,'TileLength',32);
    setTag(t,'TileWidth',32);
    setTag(t,'PlanarConfiguration',Tiff.PlanarConfiguration.Chunky);
    write(t,data);  
    close(t);