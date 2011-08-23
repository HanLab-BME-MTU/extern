function [ output ] = DVImgGetDataTypeMax( Stream )

if Stream == 0
    helpdlg('Returns the maximum value of the pixel data type the image is encoded with.  For example, an image encoded with data type 1 (2 byte, signed integer) would return a value of 32,767 (2^15-1), while an image encoded with data type 6 (2 byte, unsigned integer) would return a value of 65,535 (2^16-1).','double DVImgGetDataTypeMax(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetDataTypeMax',Stream);
end
