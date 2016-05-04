function [ output ] = DVImgGetOriginZ( Stream )

if Stream == 0
    helpdlg('Returns the Z position of the image origin in microns.','double DVImgGetOriginZ(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetOriginZ',Stream);
end
