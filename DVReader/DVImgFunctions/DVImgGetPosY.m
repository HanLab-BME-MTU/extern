function [ output ] = DVImgGetPosY( Stream,Z,W,T )

if Stream == 0
    helpdlg('Returns the Y coordinate of the section at (Z, W, T) in microns.  When the position is unknown, this function will return 0.','double DVImgGetPosY(int Stream,int Z,int W,int T)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgGetPosY',Stream,Z,W,T);
end
