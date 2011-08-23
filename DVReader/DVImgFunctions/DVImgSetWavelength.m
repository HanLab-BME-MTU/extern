function [ output ] = DVImgSetWavelength( Stream,WaveNum,Wavelength )

if Stream == 0
    helpdlg('Sets the wavelength of channel WaveNum (0 to NumW-1) within the selected image. Wavelengths are represented in nanometers.','int DVImgSetWavelength(int Stream,int WaveNum,double Wavelength)');
else
    DVImgCheckZWT(Stream,0,WaveNum,0);
    output = calllib(DVImgLibName,'DVImgSetWavelength',Stream,WaveNum,Wavelength);
    DVImgPrintErrText(output);
end
