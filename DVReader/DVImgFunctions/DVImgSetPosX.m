function [ output ] = DVImgSetPosX( Stream,Z,W,T,X )

if Stream == 0
    helpdlg('Sets the X location of the section at (Z, W, T) in microns.','int DVImgSetPosX(int Stream,int Z,int W,int T,double PosX)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgSetPosX',Stream,Z,W,T,X);
    DVImgPrintErrText(output);
end
