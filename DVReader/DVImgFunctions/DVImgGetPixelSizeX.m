function [ output ] = DVImgGetPixelSizeX( Stream )

if Stream == 0
    helpdlg('Returns the pixel size of the image in the X dimension in microns.','double DVImgGetPixelSizeX(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetPixelSizeX',Stream);
end
