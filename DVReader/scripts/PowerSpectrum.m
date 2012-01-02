function [ ] = PowerSpectrum( FileName )
%Generages a power spectrum from an image file
%This function is designed to take a file path as an input argument and
%generate the approximate power spectrum of that file and display it as a
%plot in MATLAB, eg. PowerSpectrum('/usr/local/DVImgAccess/Images/test1.dv')

% Load DV Image IO library and display version.
DVImgLibOpen(0);
libver = DVImgGetBaseLibVersion;
libbuild = DVImgGetBaseLibBuild;
fprintf('\nDVImg Library Version: %.3f, build %d\n',libver,libbuild);

DVImgOpen(1,FileName,'ro');%open the file
%get its dimensions
nx = DVImgGetNumCols(1);
ny = DVImgGetNumRows(1);
nz = DVImgGetNumZ(1);
nw = DVImgGetNumW(1);
nt = DVImgGetNumT(1);
nnx = round(nx/2);%we'll only use the first 1/4 of the image; this is approximate
nny = round(ny/2);

ReadArray = zeros(nnx,nny,nz*nw*nt);%create a 3D array to read each section into

for z=0:nz-1
    for w=0:nw-1
        for t=0:nt-1
            Section=DVImgRead(1,z,w,t); % Read section into array
            
            NewSection = fft2(Section);%take the definite fourier transform of the section
            NewSection = abs(NewSection);
            ReadArray(1:nnx,1:nny,(z+1)*(w+1)*(t+1)) = NewSection(1:nnx,1:nny);%read 1/4 of the array into ReadArray         
        end
    end
end

MaxVal = round(sqrt(nnx^2+nny^2));%this is the furthest a pixel can be from the origin
IntenTotal = zeros(2,MaxVal);%create a vector to read the intensities into
for x=1:nnx
    for y=1:nny
        Index = round(sqrt(x^2+y^2));%this is about how far a specific pixel is from the origin
        Value = mean(mean((ReadArray(x,y,1:(z+1)*(w+1)*(t+1)))));%This is the mean intensity through all the sections at that pixel
        IntenTotal(1,Index) = IntenTotal(1,Index)+Value;
        IntenTotal(2,Index) = IntenTotal(2,Index)+1;%counts how many values have been added to average them later
    end
end

IntenTotalNew = 0;
NewPlace = 0;

for place=1:MaxVal
        NewPlace = NewPlace+1;
        IntenTotalNew(NewPlace) = IntenTotal(1,place)/IntenTotal(2,place);
end

IntenTotalNew = log(IntenTotalNew);%scale it logarithmically for more convenient display
plot(IntenTotalNew);%finally, plot the data

DVImgClose(1);
DVImgLibClose();
fprintf('Done!\n');