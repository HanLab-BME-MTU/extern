function [ output ] = DVImgSetLensID( Stream,LensID )

if Stream == 0
    helpdlg('Sets the ID value of the lens that was used to acquire the image.','int DVImgSetLensID(int Stream,int LensID)');
else
    output = calllib(DVImgLibName,'DVImgSetLensID',Stream,LensID);
    DVImgPrintErrText(output);
end
