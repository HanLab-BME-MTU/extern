function [ output ] = DVImgGetOriginY( Stream )

if Stream == 0
    helpdlg('Returns the Y position of the image origin in microns.','double DVImgGetOriginY(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetOriginY',Stream);
end
