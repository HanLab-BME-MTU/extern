function [ output ] = DVImgSetOrigin( Stream,X,Y,Z )

if Stream == 0
    helpdlg('Sets the origin point of the image in microns.','int DVImgSetOrigin(int Stream,double XOrig,double YOrig,double ZOrig)');
else
    output = calllib(DVImgLibName,'DVImgSetOrigin',Stream,X,Y,Z);
    DVImgPrintErrText(output);
end
