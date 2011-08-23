function [ output ] = DVImgSetDisplay( Stream,WaveNum,Min,Max,Exp )

if Stream == 0
    helpdlg('Sets the display characteristics of the image.  ScaleMin and ScaleMax default to the minimum and maximum intensity levels of the image, but often these default values are suboptimal because a few outlying pixels can throw off the whole scale.  ScaleExp sets the base for an exponential scaling factor (for a linear scale, set it to 1).','int DVImgSetDisplay(int Stream,int WaveNum,double ScaleMin, double ScaleMax, double ScaleExp)');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgSetDisplay',Stream,WaveNum,Min,Max,Exp);
    DVImgPrintErrText(output);
end
