function [ output ] = DVImgCheckFile( FileName )

output = calllib(DVImgLibName,'DVImgCheckFile',FileName);
DVImgPrintErrText(output);
