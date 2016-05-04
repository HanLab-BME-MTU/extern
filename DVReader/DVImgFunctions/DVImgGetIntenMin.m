function [ output ] = DVImgGetIntenMin( Stream,WaveNum )

if Stream == 0
    helpdlg('Returns the minimum intensity of a specific channel (0 to NumW-1).  Note that if operations have been performed on this channel since the minimum intensity was last calculated, it may be incorrect.','double DVImgGetIntenMin(int Stream,int WaveNum)');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgGetIntenMin',Stream,WaveNum);
end
