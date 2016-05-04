function [ output ] = DVImgGetDisplayMin( Stream,WaveNum )

if Stream == 0
    helpdlg('Returns the current minimum display level for the specified image and wavelength.  If an explicit one has not been set, this defaults to the minimum intensity in that wavelength.  Any pixels at a lower intensity level than the one set here will be black in the display window.','double DVImgGetDisplayMin(int Stream, int WaveNum)');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgGetDisplayMin',Stream,WaveNum);
end
