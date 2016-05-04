function [ output ] = DVImgGetPixelSizeY( Stream )

if Stream == 0
    helpdlg('Returns the pixel size of the image in the Y dimension in microns.','double DVImgGetPixelSizeY(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetPixelSizeY',Stream);
end
