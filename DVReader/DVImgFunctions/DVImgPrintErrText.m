function [ ] = DVImgPrintErrText( ErrorNum )

if ErrorNum ~= 0
    error(DVImgGetErrText(ErrorNum));
end
