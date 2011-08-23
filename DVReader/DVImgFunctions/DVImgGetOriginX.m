function [ output ] = DVImgGetOriginX( Stream )

if Stream == 0
    helpdlg('Returns the X position of the image origin in microns.','double DVImgGetOriginX(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetOriginX',Stream);
end
