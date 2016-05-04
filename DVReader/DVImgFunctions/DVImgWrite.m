function [ output ] = DVImgWrite( Stream,Z,W,T,Array )

if Stream == 0
    helpdlg('Writes Array to the section at (Z, W, T) in the requested Stream.  Array should be an array of single precision floats with XY dimensions corresponding to those of the requested Stream.  The array will be converted to the appropriate data type while it is written to the storage device.','int DVImgWrite(int Stream, int Z, int W, int T, float Array)');
else
    Err = calllib(DVImgLibName,'DVImgMoveToZWT',Stream,Z,W,T);
    DVImgPrintErrText(Err);
    output = calllib(DVImgLibName,'DVImgWrite',Stream,Array);
    DVImgPrintErrText(output);
end
