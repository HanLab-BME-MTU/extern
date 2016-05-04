function [ output ] = DVImgMoveToZWT( Stream,Z,W,T )

if Stream == 0
    helpdlg('Moves the stream pointer to a particular (Z, W, T) coordinate.','int DVImgMoveToZWT(int Stream,int Z,int W,int T)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgMoveToZWT',Stream,Z,W,T);
    DVImgPrintErrText(output);
end
