function [ output ] = DVImgGetLensID( Stream )

if Stream == 0
    helpdlg('Returns the ID of the lens used to acquire the image. See the lens data base (dvlenses.tab) for further information.','int DVImgGetLensID(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetLensID',Stream);
end
