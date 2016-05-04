%This is a simple script that determines the correct min, mean and max
%values for a set of files listed in a text document and sets them.

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
fprintf('Testing file %d/%d\n',file,d(1));
FileName=char(FileNameCell{1}(file));%creates a character array for the file name
InFileName=strrep(FileName,' ','');%removes spaces from a line
OutFileName=strrep(InFileName,'.dv','_NEW.dv');%adds _NEW to given filename
DVImgOpen(1,InFileName,'rw'); %open existing DV image

%get the image's dimensions
nx = DVImgGetNumCols(1);
ny = DVImgGetNumRows(1);
nz = DVImgGetNumZ(1);
nw = DVImgGetNumW(1);
nt = DVImgGetNumT(1);

TypeMax = DVImgGetDataTypeMax(1);
TypeMin = DVImgGetDataTypeMin(1);

for w=0:nw-1 %wavelengths are 0-indexed
    %get the current min/max stats
    minimum = TypeMax;%the min can't be greater than the file max, so we'll start it there
    maximum = TypeMin;%the max can't be less than the file min, so we'll start it there
    total = 0;%we'll count up the means here so we can average them
    
   
    for z=0:nz-1
        for t=0:nt-1
            Section=DVImgRead(1,z,w,t); % Read section
            
            MinSection = min(min(Section));
            MaxSection = max(max(Section));
            MeanSection = mean(mean(Section));
            
            DVImgSetMin(1,z,w,t,MinSection);%set the section minimum
            DVImgSetMax(1,z,w,t,MaxSection);%set the section maximum
            DVImgSetMean(1,z,w,t,MeanSection);%set the section mean
            
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
    DVImgSetIntenStats(1,w,minimum,maximum,meanT);
    fprintf('Wavelength %d: min: %d max: %d mean: %.3f\n',w,minimum,maximum,meanT);

end
   
DVImgClose(1);

end
DVImgLibClose();

fprintf('Done!\n');
