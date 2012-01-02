%This script sets the display statistics for an image based on cutting a
%certain fraction of them off the top and the bottom.  It reads files from
%a text document called DisplayStats.txt and sets their display statistics
%according to those parameters.  This is a very simple script and probably
%does not work with floating-point data types

MinCutOff = 1/100;%set the fraction of pixels to ignore at the bottom
MaxCutOff = 1/100;%Set the fraction of pixels to ignore at the top

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
fprintf('Setting file %d/%d\n',file,d(1));
FileName=char(FileNameCell{1}(file));%creates a character array for the file name
InFileName=strrep(FileName,' ','');%removes spaces from a line
Error = DVImgCheckFile(InFileName);
if Error ~= 0 %If there is an error
    %skip this file and move on
else%continue

DVImgOpen(1,InFileName,'rw'); %open existing DV image

%get the image's dimensions and data type info
nx = DVImgGetNumCols(1);
ny = DVImgGetNumRows(1);
nz = DVImgGetNumZ(1);
nw = DVImgGetNumW(1);
nt = DVImgGetNumT(1);
TypeMax = DVImgGetDataTypeMax(1);
TypeMin = DVImgGetDataTypeMin(1);

%load each pixel's intensity into the histogram.  Also check min/max.
for w=0:nw-1
    
    Min = TypeMax;%find the min value to start polling the histogram later on
    Max = TypeMin;%find the max value to start polling the histogram later on
    Histogram = zeros(1,TypeMax - TypeMin + 1);%create a histogram we'll load the pixel data into
    
    for z=0:nz-1
        for t=0:nt-1
            Section=DVImgRead(1,z,w,t); % Read section into array
            MinSection = min(min(Section));
            MaxSection = max(max(Section));
            if MinSection < Min
                Min = min(min(Section));
            end
            if MaxSection > Max
                Max = max(max(Section));
            end
            for x=1:nx
                for y=1:ny
                    PixValue = Section(x,y);
                    Histogram(PixValue-TypeMin+1) = Histogram(PixValue-TypeMin+1)+1;
                end
            end
        end
    end
    
    MinCutOffNum = MinCutOff*nz*nt*nx*ny;%set the number of pixels to cut off the bottom
    MaxCutOffNum = MaxCutOff*nz*nt*nx*ny;%set the number of pixels to cut off the top
    MinSum = 0;
    MaxSum = 0;
    i = Min+1;
    j = Max-TypeMin+1;

    while (MinSum < MinCutOffNum)%find what pixel intensity to chop the top off at
        MinSum = Histogram(i) + MinSum;
        i = i + 1;
    end
    while (MaxSum < MaxCutOffNum)%find what pixel intensity to chop the bottom off at
        MaxSum = Histogram(j) + MaxSum;
        j = j - 1;
    end

    i = i-TypeMin;
    j = j-TypeMin;

    DVImgSetDisplay(1,w,i,j,1);
    
    Mean = DVImgGetIntenMean(1,w);%we didn't calculate mean, so we'll use the old value
    DVImgSetIntenStats(1,w,Min,Max,Mean);%set the intensity statistics, as long as we have them might as well

    fprintf('Wavelength: %d Min: %d Max: %d\n',w,i,j);
end

DVImgClose(1);
end %end conditional checking if file is error-free
end %end file for loop
DVImgLibClose();
fprintf('Done!\n');
