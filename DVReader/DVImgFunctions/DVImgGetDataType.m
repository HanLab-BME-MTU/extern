function [ output ] = DVImgGetDataType( Stream )

if Stream == 0
    helpdlg('Returns the pixel data type, 0-8.  Below is a list of each possibility and its significance.  0=integer, one byte; 1=integer, two byte signed integer; 2=floating-point, four byte; 6=unsigned integer, two byte unsigned integer; 7=integer, four byte signed integer; 8=long float, eight byte floating-point.','int DVImgGetDataType(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetDataType',Stream);
end
