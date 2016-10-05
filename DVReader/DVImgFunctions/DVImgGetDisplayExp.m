function [ output ] = DVImgGetDisplayExp( Stream,WaveNum )

if Stream == 0
    helpdlg('Returns the exponential scaling factor for displaying the specified image and wavelength.  If an explicit one has not been set, this defaults to 1 (linear).','double DVImgGetDisplayExp(int Stream, int WaveNum) ');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgGetDisplayExp',Stream,WaveNum);
end
