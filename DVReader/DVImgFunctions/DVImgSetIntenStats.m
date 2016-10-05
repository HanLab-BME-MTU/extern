function [ output ] = DVImgSetIntenStats( Stream,WaveNum,Min,Max,Mean )

if Stream == 0
    helpdlg('This sets the intensity statistics for the requested channel of the image (minimum intensity, maximum intensity and mean intensity respectively).  Note that the library will not check the intensity values.','int DVImgSetIntenStats(int Stream,int WaveNum,double Min,double Max,double Mean)');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgSetIntenStats',Stream,WaveNum,Min,Max,Mean);
    DVImgPrintErrText(output);
end
