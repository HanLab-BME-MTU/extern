function [ output ] = DVImgSetPixelSize( Stream,X,Y,Z )

if Stream == 0
    helpdlg('Sets the X, Y and Z pixel sizes in microns.','int DVImgSetPixelSize(int Stream, double DX, double DY, double DZ) ');
else
    output = calllib(DVImgLibName,'DVImgSetPixelSize',Stream,X,Y,Z);
    DVImgPrintErrText(output);
end
