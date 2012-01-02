% This program accepts a text file with the full path for 3+ channel images, subtracts background
% autofluorescence from a photoactivated image, then creates an output image with the difference channel as well as
% channels 4+ of the original image.  The program assumes the first 3
% channels of the input image are used for 1. collecting the background
% autofluorescence image, 2. photoactivating, and 3. collecting the photoactivated image.
% The text file with a list of files for batch processing should have each file name on a separate line with no extra
% spaces.  Ouput files are put in the same directory as the input files.

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
FileNameCell=textscan(fid,'%s', 'delimiter', '\n','MultipleDelimsAsOne',1);%creates a cell array of strings for the file names
d=size(FileNameCell{1});%determines size of cell array = # of files to be processed
fclose(fid); % closes text file

%Performs calculations on each file in turn.
for file=1:d(1)
fprintf('Processing file %d/%d\n',file,d(1));
FileName=char(FileNameCell{1}(file));%creates a character array for the file name
InFileName=strrep(FileName,' ','');%removes spaces from a line
OutFileName=strrep(InFileName,'.dv','_SUB.dv');%adds _SUB to given filename for the output

DVImgOpen(1,InFileName,'ro');% open the file to be processed (3+ channel DV image)

% get dimensions and pixel sizes of input image
nx = DVImgGetNumCols(1);
ny = DVImgGetNumRows(1);
nz = DVImgGetNumZ(1);
nw = DVImgGetNumW(1);
nt = DVImgGetNumT(1);
pixX = DVImgGetPixelSizeX(1);
pixY = DVImgGetPixelSizeY(1);
pixZ = DVImgGetPixelSizeZ(1);
DataType = DVImgGetDataType(1);
ImgSequence = DVImgGetImageSequence(1);
 
% Create output image file with same dimensions minus 2 wavelengths
DVImgCreate(2,OutFileName,nx,ny,nz,nw-2,nt,pixX,pixY,pixZ,ImgSequence,DataType,1);

%Get min/max/Wavelength values for each output channel
Max=zeros(nw-2);%create array for max values.  Initialized to zero.
Min=zeros(nw-2);%create array for min values.
Mean=zeros(nw-2);%create array for mean values.
Min(1)=DVImgGetIntenMax(1,2); %Initialize diff channel min to activated channel max.
Wave=zeros(nw-2);%create array for wavelength names
Wave(1)=DVImgGetWavelength(1,0);%define difference channel to have same wavelength as Background image

for w=3:nw-1
   Max(w-1)=DVImgGetIntenMax(1,w);
   Min(w-1)=DVImgGetIntenMin(1,w);
   Mean(w-1)=DVImgGetIntenMean(1,w);
   Wave(w-1)=DVImgGetWavelength(1,w);
end

% 1 section at a time, subtract background autofluorescence from activated channel 
     for t=0:nt-1
         for z=0:nz-1
                 Section1=DVImgRead(1,z,0,t); % Read section
                 Section2=DVImgRead(1,z,2,t); % Read section
                 Section3=Section2-Section1;%subtract background level from photoactivated level
                 DVImgWrite(2,z,0,t,Section3); % Write difference channel in output image
         end
     end

    
     %Copy any additional channels into output image
     for t=0:nt-1
         for w=3:nw-1
             for z=0:nz-1
                 Section1=DVImgRead(1,z,w,t); % Read section
                 DVImgWrite(2,z,w-2,t,Section1); % copy channel to output image
             end
         end
     end

 %Write Min/Max/wavelength values for each channel into output image
 for m=0:nw-3
    DVImgSetIntenStats(2,m,Min(m+1),Max(m+1),Mean(m+1));
    DVImgSetWavelength(2,m,Wave(m+1));
 end 

%Write display information for channels 4+ to channels 2+
for w=3:nw-1
    DisplayMin = DVImgGetDisplayMin(1,w);
    DisplayMax = DVImgGetDisplayMax(1,w);
    DisplayExp = DVImgGetDisplayExp(1,w);
    DVImgSetDisplay(2,w-2,DisplayMin,DisplayMax,DisplayExp);
end
 
%copy the extended header for channels 3+ to 1+
for z=0:nz-1
    for w=2:nw-1
        for t=0:nt-1
            PosX = DVImgGetPosX(1,z,w,t);
            PosY = DVImgGetPosY(1,z,w,t);
            PosZ = DVImgGetPosZ(1,z,w,t);
            Time = DVImgGetTime(1,z,w,t);
            PhotoVal = DVImgGetPhotoVal(1,z,w,t);
            Min = DVImgGetMin(1,z,w,t);
            Max = DVImgGetMax(1,z,w,t);
            Mean = DVImgGetMean(1,z,w,t);
             
            DVImgSetPosX(2,z,w-2,t,PosX);
            DVImgSetPosY(2,z,w-2,t,PosY);
            DVImgSetPosZ(2,z,w-2,t,PosZ);
            DVImgSetTime(2,z,w-2,t,Time);
            DVImgSetPhotoSen(2,z,w-2,t,PhotoVal);
            DVImgSetMin(2,z,w-2,t,Min);
            DVImgSetMax(2,z,w-2,t,Max);
            DVImgSetMean(2,z,w-2,t,Mean);
        end
    end
end


%Channel 0 will need its min/max/mean recalculated for each sec and total
minimum = DVImgGetDataTypeMax(2);
maximum = DVImgGetDataTypeMin(2);

total = 0;

for z=0:nz-1
    for t=0:nt-1
        Section1=DVImgRead(2,z,0,t); % Read section
        
        MinSection = min(min(Section1));
        MaxSection = max(max(Section1));
        MeanSection = mean(mean(Section1));
       
        DVImgSetMin(1,z,0,t,MinSection);%set the section minimum
        DVImgSetMax(1,z,0,t,MaxSection);%set the section maximum
        DVImgSetMean(1,z,0,t,MeanSection);%set the section mean
        
        total = total + MeanSection;
        if MinSection < minimum
            minimum = MinSection;
        end
        if MaxSection > maximum
            maximum = MaxSection;
        end
    end
end

 
DVImgClose(1); %Close input image
DVImgClose(2); % Close output image
end %for loop for file names
DVImgLibClose();
fprintf('Done \n');
