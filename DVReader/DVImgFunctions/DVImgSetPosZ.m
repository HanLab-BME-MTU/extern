function [ output ] = DVImgSetPosZ( Stream,Z,W,T,PosZ )

if Stream == 0
    helpdlg('Sets the Z location of the section at (Z, W, T) in microns.','int DVImgSetPosZ(int Stream,int Z,int W,int T,double PosZ)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgSetPosZ',Stream,Z,W,T,PosZ);
    DVImgPrintErrText(output);
end
