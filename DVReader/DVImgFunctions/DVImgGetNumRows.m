function [ output ] = DVImgGetNumRows( Stream )

if Stream == 0
    helpdlg('Returns the number of rows (y dimension) in an image.','int DVImgGetNumRows(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetNumRows',Stream);
end
