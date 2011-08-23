function [ output ] = DVImgSetMean( Stream,Z,W,T,Mean )

if Stream == 0
    helpdlg('Sets the mean intensity of the section at (Z, W, T).  This is useful if an operation has been performed on this section that might have altered the mean intensity.','int DVImgSetMean(int Stream,int Z,int W,int T,double Mean)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgSetMean',Stream,Z,W,T,Mean);
    DVImgPrintErrText(output);
end
