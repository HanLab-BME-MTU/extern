%This takes a list of input images and transforms them all to either ZTW,
%WZT, or ZWT.  Although most software will read the images in any of the
%above orders, it can be advantageous to format them in the order they will
%be read for faster loading from the hard drive.

Order = 'ZTW'; %Put ZTW, WZT, or ZWT here depending on how you want your output images formatted

%Get the number the order corresponds to for creating the file
if strcmpi('ZTW',Order)
    OrderNum = 0;
elseif strcmpi('WZT',Order)
    OrderNum = 1;
elseif strcmpi('ZWT',Order)
    OrderNum = 2;
end

% Load DV Image IO library and display version.
DVImgLibOpen(0);
libver = DVImgGetBaseLibVersion;
libbuild = DVImgGetBaseLibBuild;
fprintf('\nDVImg Library Version: %.3f, build %d\n',libver,libbuild);

%change to image directory and read text file with file names to be processed
cd (dvimgroot); %Change directory to dvimgroot directory
cd 'Images'; %Change directory to Sample Images directory
batchfile=fullfile(dvimgroot, 'Scripts', 'batch_process.txt');
fid=fopen(batchfile);%opens text file with file names to be processed
FileNameCell=textscan(fid,'%s', 'delimiter', '\n','MultipleDelimsAsOne',1);%creates a cell array of strings
d=size(FileNameCell{1});%determines size of cell array = # of files to be processed
fclose(fid);

%Performs calculations on each file in turn.
for file=1:d(1)
fprintf('Converting file %d/%d\n',file,d(1));
FileName=char(FileNameCell{1}(file));%creates a character array for the file name
InFileName=strrep(FileName,' ','');%removes spaces from a line
Error = DVImgCheckFile(InFileName);
if Error ~= 0 %If there is an error
    %skip this file and move on
else%continue
OutFileName=strrep(InFileName,'.dv','_TMPO.dv');%adds _TMPO (temporary) to given filename
OutFileName=strrep(OutFileName,'TMPO',Order);%adds the order suffix
DVImgOpen(1,InFileName,'ro'); %open existing DV image

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
ChannelArray = zeros(nw,7);
for Channel=0:nw-1
    index = Channel+1;
    Wavelength = DVImgGetWavelength(1,Channel);
    IntenMin = DVImgGetIntenMin(1,Channel);
    IntenMax = DVImgGetIntenMax(1,Channel);
    IntenMean = DVImgGetIntenMean(1,Channel);
    ChannelArray(index,1) = Wavelength;
    ChannelArray(index,2) = IntenMin;
    ChannelArray(index,3) = IntenMax;
    ChannelArray(index,4) = IntenMean;
    ChannelArray(index,5) = DVImgGetDisplayMin(1,Channel);
    ChannelArray(index,6) = DVImgGetDisplayMax(1,Channel);
    ChannelArray(index,7) = DVImgGetDisplayExp(1,Channel);
end
    
%create and sync the header
DVImgCreate(2,OutFileName,nx,ny,nz,nw,nt,pixX,pixY,pixZ,OrderNum,DataType,1);
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

%sync the data
for t=0:nt-1
    for w=0:nw-1
        for z=0:nz-1
            Section=DVImgRead(1,z,w,t); % Read section into array
            DVImgWrite(2,z,w,t,Section); % copy channel to output image
        end
    end
end

DVImgClose(1);
DVImgClose(2);
end %end conditional checking if file is error-free
end %end file for loop

DVImgLibClose;
fprintf('Done!\n');
