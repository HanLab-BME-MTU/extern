function [ output ] = DVImgGetNumW( Stream )

if Stream == 0
    helpdlg('Returns the number of channels (wavelengths) in an image.','int DVImgGetNumW(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetNumW',Stream);
end
