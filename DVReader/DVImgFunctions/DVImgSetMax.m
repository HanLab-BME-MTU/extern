function [ output ] = DVImgSetMax( Stream,Z,W,T,Max )

if Stream == 0
    helpdlg('Sets the maximum intensity of the section at (Z, W, T).  This is useful if an operation has been performed on this section that might have altered the maximum value.','int DVImgSetMax(int Stream,int Z,int W,int T,double Max)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgSetMax',Stream,Z,W,T,Max);
    DVImgPrintErrText(output);
end
