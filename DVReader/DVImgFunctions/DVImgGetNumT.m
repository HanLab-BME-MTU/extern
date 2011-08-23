function [ output ] = DVImgGetNumT( Stream )

if Stream == 0
    helpdlg('Returns the number of time points in an image.','int DVImgGetNumT(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetNumT',Stream);
end
