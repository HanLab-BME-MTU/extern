function [ output ] = DVImgGetLibVersion( )

output = calllib(DVImgLibName,'DVImgGetLibVersion');
