function [ output ] = DVImgGetPixelSizeZ( Stream )

if Stream == 0
    helpdlg('Returns the pixel size of the image in the Z dimension in microns.','double DVImgGetPixelSizeZ(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetPixelSizeZ',Stream);
end
