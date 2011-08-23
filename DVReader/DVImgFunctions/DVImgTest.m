function [ output ] = DVImgTest( ImgPath )

% Load DV Image IO library and display version.
libver = DVImgGetBaseLibVersion;
libbuild = DVImgGetBaseLibBuild;
fprintf('\nDVImg Library Version: %.3f, build %d\n',libver,libbuild);

%Check to be sure the file exists
DVImgCheckFile(ImgPath);

%loads an image
DVImgOpen(1,ImgPath,'rw');
    
% Gets everything we know about the pixels
DataType = DVImgGetDataType(1);
Bytes = DVImgGetBytesPerPixel(1);
MinValue = DVImgGetDataTypeMin(1);
MaxValue = DVImgGetDataTypeMax(1);
fprintf('\nPixel Data Type: %d\nBytes per pixel: %d\nMinimum possible intensity: %d\nMaximum possible intensity: %d\n',DataType,Bytes,MinValue,MaxValue);

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
fprintf('Writable range: 0-%d(Z), 0-%d(W), 0-%d(T)\n',NumZ-1,NumW-1,NumT-1);

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
    fprintf('Channel %d: w/l: %dnm  INTENSITY: min:%d max:%d mean:%.3f  DISPLAY: min:%d max:%d exp:%d\n',Channel,Wavelength,IntenMin,IntenMax,IntenMean,DispMin,DispMax,DispExp);
end