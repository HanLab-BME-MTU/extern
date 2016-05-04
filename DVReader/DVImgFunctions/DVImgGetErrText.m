function [ output ] = DVImgGetErrText( ErrorNum )

output = calllib(DVImgLibName,'DVImgGetErrText',ErrorNum);
