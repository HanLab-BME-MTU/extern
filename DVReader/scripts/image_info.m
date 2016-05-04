% This program reads an image and outputs diagnostic information about it.
% It is a practical demonstration of many of the "Get" commands.

% Load DV Image IO library and display version.
DVImgLibOpen(0);
libver = DVImgGetBaseLibVersion;
libbuild = DVImgGetBaseLibBuild;
fprintf('\nDVImg Library Version: %.3f, build %d\n',libver,libbuild);
fprintf('\nPlease select the image you would like to analyze\n');
[fname, fpath] = uigetfile; %enter the full path of the image you want to analyze at the prompt

%loads an image
DVImgOpen(1,[fpath fname],'ro');
    
% Gets everything we know about the pixels
DataType = DVImgGetDataType(1);
Bytes = DVImgGetBytesPerPixel(1);
MinValue = DVImgGetDataTypeMin(1);
MaxValue = DVImgGetDataTypeMax(1);
fprintf('\nPixel Data Type: %d\nBytes per pixel: %d\nMinimum possible intensity: %d\nMaximum possible intensity: %d\n',DataType,Bytes,MinValue,MaxValue);

%Gets the amount of space in the extended header
HeaderSpace = DVImgGetFieldHdrSpace(1);
fprintf('\nFree space in the extended header: %d\n',HeaderSpace);

%Get the image sequence
ImgSequence = DVImgGetImageSequence(1);
if ImgSequence == 0
    ImgString = 'ZTW';
elseif ImgSequence == 1
    ImgString = 'WZT';
elseif ImgSequence == 2
    ImgString = 'ZWT';
end
DispString = ['\nSection sequence: ',ImgString,'\n'];
fprintf(DispString);

% Gets the XYZ origins
OriginX = DVImgGetOriginX(1);
OriginY = DVImgGetOriginY(1);
OriginZ = DVImgGetOriginZ(1);
fprintf('\nXYZ origin point (in um): (%d, %d, %d)\n',OriginX,OriginY,OriginZ);

% Gets the Lens ID the image was created with
LensID = DVImgGetLensID(1);
fprintf('\nLens ID: %d\n', LensID);

% Gets the image's dimensions
NumX = DVImgGetNumCols(1);
NumY = DVImgGetNumRows(1);
NumZ = DVImgGetNumZ(1);
NumW = DVImgGetNumW(1);
NumT = DVImgGetNumT(1);
fprintf('\nDimensions (pixels):  %d(X) x %d(Y) x %d(Z) x %d(Wavelengths) x %d(Time points)\n',NumX,NumY,NumZ,NumW,NumT);

% Gets the pixel sizes and calculates the image's real dimensions
ScaleX = DVImgGetPixelSizeX(1);
ScaleY = DVImgGetPixelSizeY(1);
ScaleZ = DVImgGetPixelSizeZ(1);
LengthX = NumX*ScaleX;
LengthY = NumY*ScaleY;
LengthZ = NumZ*ScaleZ;
fprintf('Dimensions (microns): %.3f(X) x %.3f(Y) x %.3f(Z)\n',LengthX,LengthY,LengthZ);

% Gets the min/max/mean image intensities and wavelength for each channel
fprintf('\n');
for Channel = 0:NumW-1
    Wavelength = DVImgGetWavelength(1,Channel);
    IntenMin = DVImgGetIntenMin(1,Channel);
    IntenMax = DVImgGetIntenMax(1,Channel);
    IntenMean = DVImgGetIntenMean(1,Channel);
    DispMin = DVImgGetDisplayMin(1,Channel);
    DispMax = DVImgGetDisplayMax(1,Channel);
    DispExp = DVImgGetDisplayExp(1,Channel);
    fprintf('Channel %d: Wavelength %dnm\n  MinInten: %d  MaxInten: %d  MeanInten: %.3f  MinDisp: %d  MaxDisp: %d  ExpDisp: %d\n',Channel,Wavelength,IntenMin,IntenMax,IntenMean,DispMin,DispMax,DispExp);
end
DVImgLibClose();