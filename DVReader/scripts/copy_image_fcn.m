%This function does the same thing as the copy_image.m file, except it can
%be called as a function from elsewhere.  To use it, simply type
%copy_img_fcn, and in parenthesis following it, pass the full path name of the
%file you wish to copy.

function [] = copy_image_fcn(FileName)

DVImgLibOpen(0);
libver = DVImgGetBaseLibVersion;
libbuild = DVImgGetBaseLibBuild;
fprintf('\nDVImg Library Version: %.3f, build %d\n',libver,libbuild);

DVImgOpen(1,FileName,'ro'); %open existing DV image
OutFileName=strrep(FileName,'.dv','_COPY.dv');%adds _COPY to given filename

% Read some parameters we'll use when creating and syncing the copy
nx = DVImgGetNumCols(1);
ny = DVImgGetNumRows(1);
nz = DVImgGetNumZ(1);
nw = DVImgGetNumW(1);
nt = DVImgGetNumT(1);
pixX = DVImgGetPixelSizeX(1);
pixY = DVImgGetPixelSizeY(1);
pixZ = DVImgGetPixelSizeZ(1);
DataType = DVImgGetDataType(1);
OriginX = DVImgGetOriginX(1);
OriginY = DVImgGetOriginY(1);
OriginZ = DVImgGetOriginZ(1);
LensID = DVImgGetLensID(1);
TimeStamp = DVImgGetTimeStamp(1);
ImgSequence = DVImgGetImageSequence(1);
ChannelArray = zeros(nw,7);
for Channel=0:nw-1
    index = Channel+1;
    ChannelArray(index,1) = DVImgGetWavelength(1,Channel);
    ChannelArray(index,2) = DVImgGetIntenMin(1,Channel);
    ChannelArray(index,3) = DVImgGetIntenMax(1,Channel);
    ChannelArray(index,4) = DVImgGetIntenMean(1,Channel);
    ChannelArray(index,5) = DVImgGetDisplayMin(1,Channel);
    ChannelArray(index,6) = DVImgGetDisplayMax(1,Channel);
    ChannelArray(index,7) = DVImgGetDisplayExp(1,Channel);
end
    
%create and sync the header
DVImgCreate(2,OutFileName,nx,ny,nz,nw,nt,pixX,pixY,pixZ,ImgSequence,DataType,1);
DVImgSetOrigin(2,OriginX,OriginY,OriginZ);
DVImgSetTimeStamp(2,TimeStamp);
DVImgSetLensID(2,LensID);
for Channel=0:nw-1
    index = Channel+1;
    DVImgSetWavelength(2,Channel,ChannelArray(index,1));
    DVImgSetIntenStats(2,Channel,ChannelArray(index,2),ChannelArray(index,3),ChannelArray(index,4));
    DVImgSetDisplay(2,Channel,ChannelArray(index,5),ChannelArray(index,6),ChannelArray(index,7));
end

%sync the extended header
for t=0:nt-1
    for w=0:nw-1
        for z=0:nz-1
            PosX = DVImgGetPosX(1,z,w,t);
            PosY = DVImgGetPosY(1,z,w,t);
            PosZ = DVImgGetPosZ(1,z,w,t);
            Time = DVImgGetTime(1,z,w,t);
            PhotoVal = DVImgGetPhotoVal(1,z,w,t);
            Min = DVImgGetMin(1,z,w,t);
            Max = DVImgGetMax(1,z,w,t);
            Mean = DVImgGetMean(1,z,w,t);
            
            DVImgSetPosX(2,z,w,t,PosX);
            DVImgSetPosY(2,z,w,t,PosY);
            DVImgSetPosZ(2,z,w,t,PosZ);
            DVImgSetTime(2,z,w,t,Time);
            DVImgSetPhotoSen(2,z,w,t,PhotoVal);
            DVImgSetMin(2,z,w,t,Min);
            DVImgSetMax(2,z,w,t,Max);
            DVImgSetMean(2,z,w,t,Mean);
        end
    end
end

%sync any metadata stored in the extended header
Fields = DVImgGetExtInfo(1);
FieldNames = fieldnames(Fields);
for Field=FieldNames'
    Value = DVImgGetExtHeaderField(1,char(Field));
    Type = char(Fields.(char(Field))(1));
    DVImgSetExtHeaderField(2,char(Field),Type,Value);
end

%sync the data
for t=0:nt-1
    for w=0:nw-1
        for z=0:nz-1
            Section = DVImgRead(1,z,w,t); % Read section
            DVImgWrite(2,z,w,t,Section); % copy section to output image
        end
    end
end

DVImgClose(1);
DVImgClose(2);
DVImgLibClose();
fprintf('Done!\n');


