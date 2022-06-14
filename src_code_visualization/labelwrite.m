function labelwrite(data, label, path)
    outIm = imdisplayWithMapColor4D(data,label);
    write4dTiffRGB(uint8(outIm),[path '.tif']);
end