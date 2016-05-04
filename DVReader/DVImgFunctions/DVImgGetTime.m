function [ output ] = DVImgGetTime( Stream,Z,W,T )

if Stream == 0
    helpdlg('Returns the elapsed time (in seconds) between the initial time of the entire image and the time of section (Z, W, T). ','double DVImgGetTime(int Stream,int Z,int W,int T)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgGetTime',Stream,Z,W,T);
end
