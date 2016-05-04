function [ output ] = DVImgClose( Stream )

if Stream == 0
    helpdlg('Closes the designated image stream. Closing image streams is essential. Failure to close a stream before opening a new one in the same stream or closing the library, will result in lost changes to the image header.','int DVImgClose(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgClose',Stream);
    DVImgPrintErrText(output);
end
