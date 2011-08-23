function [ output ] = DVImgRead( Stream,Z,W,T )

if Stream == 0
    helpdlg('Returns the image section at (Z, W, T) from the requested Stream. The image is always returned as an array of single precision floats, regardless of what type is contained within the stored image. That is, the library always converts images to single precision floating-point.','Array DVImgRead(int Stream,int Z,int W,int T)');
else
    Err = calllib(DVImgLibName,'DVImgMoveToZWT',Stream,Z,W,T);
    DVImgPrintErrText(Err);

    Rows = DVImgGetNumRows(Stream);
    Cols = DVImgGetNumCols(Stream);
    Section = zeros(Rows,Cols);
    pSection = libpointer('singlePtr',Section);
    calllib(DVImgLibName,'DVImgRead',Stream,pSection);
    output = get(pSection,'value');
end
