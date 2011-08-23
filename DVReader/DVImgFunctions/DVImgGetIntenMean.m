function [ output ] = DVImgGetIntenMean( Stream,WaveNum )

if Stream == 0
    helpdlg('Returns the mean intensity of a specific channel (0 to NumW-1).  Note that if operations have been performed on this channel since the mean intensity was last calculated, it may be incorrect. ','double DVImgGetIntenMean(int Stream,int WaveNum) ');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgGetIntenMean',Stream,WaveNum);
end
