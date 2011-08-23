function [ output ] = DVImgGetPosZ( Stream,Z,W,T )

if Stream == 0
    helpdlg('Returns the Z coordinate of the section at (Z, W, T) in microns. When the position is unknown, this function will return 0.','double DVImgGetPosZ(int Stream,int Z,int W,int T) ');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgGetPosZ',Stream,Z,W,T);
end
