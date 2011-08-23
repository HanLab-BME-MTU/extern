function [ output ] = DVImgGetBytesPerPixel( Stream )

if Stream == 0
    helpdlg('Returns the number of bytes per pixel in the image (1, 2, 4 or 8).','int DVImgGetBytesPerPixel(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetBytesPerPixel',Stream);
end
