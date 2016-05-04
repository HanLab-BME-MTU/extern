function [ output ] = DVImgOpen( Stream,FileName,AccessType )

if Stream == 0
    helpdlg('Loads an image to be read or manipulated.  All image functions after this point on a particular stream use the image loaded here, until DVImgClose is run.  Stream gives the file an identity that can be used in future functions to refer to it (see footnote 2).  FileName should be the full path to the file.  AccessType should be either ro (read only) or rw (read-writable).  Any changes you make to a file opened as ro will not be saved to disk when you close the image.','int DVImgOpen(int Stream, char FileName, char AccessType) ');
else
    output = calllib(DVImgLibName,'DVImgOpen',Stream,FileName,AccessType);
    DVImgPrintErrText(output);
end
