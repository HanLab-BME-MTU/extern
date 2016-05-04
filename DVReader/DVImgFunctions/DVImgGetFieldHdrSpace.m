function [ output ] = DVImgGetFieldHdrSpace( Stream )

if Stream == 0;
    helpdlg('Returns the amount of space available in the image''s extended header. Useful when determining whether there is space available to add another field.','int DVImgGetFieldHdrSpace(int iStream)');
else
    output = calllib(DVImgLibName,'DVImgGetFieldHdrSpace',Stream);
end