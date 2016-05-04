%This function runs the fourier transform on an input image and saves the
%output image in the same directory with the suffix "_FOURIER".  It uses
%the default algorithm which is slow, but works.
%F(u,v) = SUM{ f(x,y)*exp(-j*2*pi*(u*x+v*y)/N) } where you have a square
%image and N = DIM/2

File = fullfile(dvimgroot, 'Images', 'test1.dv');%change this input file as desired

% Load DV Image IO library and display version.
DVImgLibOpen(0);
libver = DVImgGetBaseLibVersion;
libbuild = DVImgGetBaseLibBuild;
fprintf('\nDVImg Library Version: %.3f, build %d\n',libver,libbuild);

Error = DVImgCheckFile(File);%check file for errors
if Error ~= 0 %If there is an error
    %skip this file and move on
else%continue

OutFileName=strrep(File,'.dv','_FOURIER.dv');%adds _FOURIER to given filename
DVImgOpen(1,File,'ro'); %open existing DV image

% Read some parameters we'll use when creating and syncing the transform
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
ChannelArray = zeros(nw);
for Channel=0:nw-1
    index = Channel+1;
    Wavelength = DVImgGetWavelength(1,Channel);
    ChannelArray(index) = Wavelength;
end

%create and sync the header
DVImgCreate(2,OutFileName,nx,ny,nz,nw,nt,pixX,pixY,pixZ,ImgSequence,DataType,1);
DVImgSetOrigin(2,OriginX,OriginY,OriginZ);
DVImgSetTimeStamp(2,TimeStamp);
DVImgSetLensID(2,LensID);
for Channel=0:nw-1
    index = Channel+1;
    DVImgSetWavelength(2,Channel,ChannelArray(index));
end

%Generate the transform, one section at a time
for t=0:nt-1
    for w=0:nw-1
        for z=0:nz-1
            Section=DVImgRead(1,z,w,t); % Read section into pointer array
            
            NewSection = fft2(Section);
            NewSection = abs(NewSection);
            %NewSection = fftshift(NewSection);%Changes the origin of the fourier transform to the center of the image
            DVImgWrite(2,z,w,t,NewSection); % copy section to output image
        end
    end
end

%regenerate the statistics for each section
TypeMax = DVImgGetDataTypeMax(2);
TypeMin = DVImgGetDataTypeMin(2);

for w=0:nw-1 %wavelengths are 0-indexed
    %get the current min/max stats
    minimum = TypeMax;%the min can't be greater than the file max, so we'll start it there
    maximum = TypeMin;%the max can't be less than the file min, so we'll start it there
    total = 0;%we'll count up the means here so we can average them
    
   
    for z=0:nz-1
        for t=0:nt-1
            Section=DVImgRead(2,z,w,t); % Read section
            
            MinSection = min(min(Section));
            MaxSection = max(max(Section));
            MeanSection = mean(mean(Section));
            
            DVImgSetMin(2,z,w,t,MinSection);%set the section minimum
            DVImgSetMax(2,z,w,t,MaxSection);%set the section maximum
            DVImgSetMean(2,z,w,t,MeanSection);%set the section mean
            
            total = total + MeanSection;
            if MinSection < minimum
                minimum = MinSection;
            end
            if MaxSection > maximum
                maximum = MaxSection;
            end
        end
    end
    
    meanT = sum(sum(total))/(nz*nt); %the mean is the average of all the section means
    DVImgSetIntenStats(2,w,minimum,maximum,meanT);

end

DVImgClose(1);
DVImgClose(2);
end %end conditional checking if file is error-free
DVImgLibClose();
fprintf('Done!\n');
