function [ output ] = DVImgSetPhotoSen( Stream,Z,W,T,PhotoSen )

if Stream == 0
    helpdlg('Sets the photosensor value for the section at (Z, W, T).  This is useful for normalizing pixel intensity values.','int DVImgSetPhotoSen(int Stream,int Z,int W,int T,double PhotoVal)');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgSetPhotoSen',Stream,Z,W,T,PhotoSen);
    DVImgPrintErrText(output);
end
