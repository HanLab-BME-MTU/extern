function [ output ] = DVImgSetMin( Stream,Z,W,T,Min )

if Stream == 0
    helpdlg('Sets the minimum intensity of the section at (Z, W, T).  This is useful if an operation has been performed on this section that may have altered the minimum intensity.','int DVImgSetMin(int Stream,int Z,int W,int T,double Min)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgSetMin',Stream,Z,W,T,Min);
    DVImgPrintErrText(output);
end
