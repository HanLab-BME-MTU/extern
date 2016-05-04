function [ output ] = DVImgSetTime( Stream,Z,W,T,dT )

if Stream == 0
    helpdlg('Sets the time offset in seconds between when the image was initiated and when the section at (Z, W, T) was scanned.','int DVImgSetTime(int Stream,int Z,int W,int T,double T)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgSetTime',Stream,Z,W,T,dT);
    DVImgPrintErrText(output);
end
