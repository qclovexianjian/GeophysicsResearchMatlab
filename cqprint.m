function cqprint( fileName )
% Print function: prints to cmyk tiff and eps files

% Print two files, one is tiff and another is eps
print(fileName,'-dtiffn','-r600');
% Convert the printed one into CMYK format
rgb = imread([fileName,'.tif']);
cform = makecform('srgb2cmyk');
lab = applycform(rgb,cform); 
imwrite(lab,[fileName,'.tif']);

% Print to eps format, still in cmyk format
print(fileName,'-deps','-cmyk','-tiff');


end

