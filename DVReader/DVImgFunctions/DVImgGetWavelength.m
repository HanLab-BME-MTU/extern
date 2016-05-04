function [ output ] = DVImgGetWavelength( Stream,WaveNum )

if Stream == 0
    helpdlg('Returns the wavelength of the specified channel (0 to NumW-1) in nanometers.','double DVImgGetWavelength(int Stream,int WaveNum)');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgGetWavelength',Stream,WaveNum);
end
