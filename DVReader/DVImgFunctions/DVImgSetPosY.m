function [ output ] = DVImgSetPosY( Stream,Z,W,T,Y )

if Stream == 0
    helpdlg('Sets the Y location of the section at (Z, W, T) in microns.','int DVImgSetPosY(int Stream,int Z,int W,int T,double PosY)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgSetPosY',Stream,Z,W,T,Y);
    DVImgPrintErrText(output);
end
