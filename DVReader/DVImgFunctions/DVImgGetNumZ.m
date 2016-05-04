function [ output ] = DVImgGetNumZ( Stream )

if Stream == 0
    helpdlg('Returns the number of Z sections of an image.','int DVImgGetNumZ(int Stream) ');
else
    output = calllib(DVImgLibName,'DVImgGetNumZ',Stream);
end
