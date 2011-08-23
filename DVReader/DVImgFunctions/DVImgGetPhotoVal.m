function [ output ] = DVImgGetPhotoVal( Stream,Z,W,T )

if Stream == 0
    helpdlg('Returns the photosensor value of the section at (Z, W, T).  This is useful for normalizing image intensities.','double DVImgGetPhotoVal(int Stream,int Z,int W,int T) ');
else
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgGetPhotoVal',Stream,Z,W,T);
end
