function [ output ] = DVImgGetDisplayMax( Stream,WaveNum )

if Stream == 0
    helpdlg('Returns the current maximum display level for the specified image and wavelength.  If an explicit one has not been set, this defaults to the maximum intensity in that wavelength.  Any pixels at a higher intensity level than the one set here will be fully saturated in the display window.','double DVImgGetDisplayMax(int Stream, int WaveNum)');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgGetDisplayMax',Stream,WaveNum);
end
