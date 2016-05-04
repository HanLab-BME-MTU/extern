function [ output ] = DVImgGetNumCols( Stream )

if Stream == 0
    helpdlg('Returns the number of columns (x dimension) in an image.','int DVImgGetNumCols(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetNumCols',Stream);
end
