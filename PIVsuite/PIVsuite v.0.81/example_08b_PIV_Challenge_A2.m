%% Example 08b - Treating test case A2 from 3rd PIV Challenge
% This case is essentially the same as example 07: The test case A2 from 3rd PIV challenge [ref. 6] is
% treated. Mean and rms velocities are computed, velocity PDF is determined and wavenumber spectra is
% calculated; these results can be compared to results shown in reference [6], Figs. 3b, 12, 14 and Table 11.
%
% *Reference:* [6] Stanislas, M., K. Okamoto, C. J. Kahler, J. Westerweel and F. Scarano, (2008): Main results
% of the third international PIV Challenge. Experiments in Fluids, vol. 45, pp. 27-71.
%
% *Instructions:*
%
% # Download images (case A2) from <http://www.pivchallenge.org/pub05/A/A2.zip PIVchallenge site>,
% # Unzip them to folder |../Data/Test PIVChallenge3A2|,
% # Run this example.

%% Define image pairs to treat
% Initialize the variable |pivPar|. Get list of image files.

% Initialize variables
clear;
pivPar = [];                                              % variable for settings
imagePath = ['..',filesep,'Data',filesep,'Test PIVChallenge3A2'];  % folder with processed images

try
    % get list of "A" images
    aux = dir([imagePath, filesep, '*a.tif']);
    for kk = 1:numel(aux), fileList{kk} = [imagePath, filesep, aux(kk).name]; end   %#ok<SAGROW>
    im1 = sort(fileList);
    
    % get list of "B" images
    aux = dir([imagePath, filesep, '*b.tif']);
    for kk = 1:numel(aux), fileList{kk} = [imagePath, filesep, aux(kk).name]; end   %#ok<SAGROW>
    im2 = sort(fileList);
catch
    % if error on loading images, probably images were not downloaded. Give a message.
    error(['No images found. Please, download images (case A2) from http://www.pivchallenge.org/pub05/A/A2.zip, '...
        'unzip them and place them to folder ../Data/Test PIVChallenge3A2.']);
end

%% Running first 4 passes of PIV processing
% First, image pairs will be treated with 4-pass PIV with decreasing size of interrogation area and spacing of
% interrogation grid. Set the corresponding parameters:

pivPar.iaSizeX = [64 32 16  8];           % interrogation area size for three passes
pivPar.iaStepX = [32 16  8  4];           % grid spacing for three passes
pivPar.anVelocityEst = 'none';            % iterate each image pair zero  
pivPar.anOnDrive = true;                  % files with results will be stored in an output folder
pivPar.anForceProcessing = false;         % if false, only image pairs, for which no file with results is 
        % available, will be processed. Processing is skipped if file with results is available. If true,
        % processing is carried out even if result file is present. (Set this parameter to true if all image
        % pairs should be reprocessed, for example because of different setting of processing parameters).
pivPar.anTargetPath = [imagePath,filesep,'pivOut1'];   % folder for storing results of first three passes
pivPar.qvPair = {'U'};
figure(1);
[pivPar1] = pivParams([],pivPar,'defaultsSeq');
            
% Run the analysis now
[pivData1] = pivAnalyzeImageSequence(im1,im2,[],pivPar1);


%% Running additional (5th) pass of PIV processing
% Images will be processed with  additional pass, during which the interrogation area size is kept at
% 8x8 pixels and grid spacing at 8x8. Results will be stored in another variable, and files will be placed
% in different folder.

pivPar.iaSizeX = [ 8];              % interrogation area size for fourth and fifth pass
pivPar.iaStepX = [ 4];              % grid spacing 
pivPar.anVelocityEst = 'pivData';   % use velocity data stored in pivData as velocity estimate used for image 
                                    % deformation. By this setting, results of previous passes are transfered
                                    % to additional passses.
pivPar.ccMethod = 'dcn';            % set cross-correlation method to direct convolution (faster than default
                                    % FFT, if displacements are small, which is the case of final passes)
pivPar.anTargetPath = [imagePath,filesep,'pivOut2']; % set output folder different than that for first 16x16 px
pivPar.qvPair = {'U'};
figure(1);
[pivPar2] = pivParams([],pivPar,'defaultsSeq');

% Run the analysis now. Note that pivData1 (results of first three passes) are an input variable.
[pivData2] = pivAnalyzeImageSequence(im1,im2,pivData1,pivPar2);

%% Show results
% First, samples of the velocity field are shown
Nt = size(pivData2.U,3);
figure(1);
for kk=1:10:Nt
    hold off;
    pivQuiver(pivData2,'timeSlice',kk,...
        'Umag');
    title('Background: U_{mag}. Quiver: velocity (black: valid, white: replaced)');
end

%%
% Crop two the image for a few pixels where some NaN's are present
pivData1 = pivManipulateData('limitX',pivData1,[10, Inf]);
pivData2 = pivManipulateData('limitX',pivData2,[10, Inf]);
N1 = numel(pivData1.U);
N2 = numel(pivData2.U);

%% 
% Compute the mean and rms velocities, display then and compare them to the reference values.

% Statistics for the velocity fields after 4 passes:
meanU1 = mean(reshape(pivData1.U,N1,1));
meanV1 = mean(reshape(pivData1.V,N1,1));
stdU1 = std(reshape(pivData1.U,N1,1));
stdV1 = std(reshape(pivData1.V,N1,1));

% Statistics for the velocity fields after 5th pass:
meanU2 = mean(reshape(pivData2.U,N2,1));
meanV2 = mean(reshape(pivData2.V,N2,1));
stdU2 = std(reshape(pivData2.U,N2,1));
stdV2 = std(reshape(pivData2.V,N2,1));

% Print results:
fprintf('\nStatistics (16x16 px): mean(U) = %+6.4f, mean(V) = %+6.4f, std(U) = %+6.4f, std(V) = %+6.4f\n',...
    meanU1,-meanV1,stdU1,stdV1);
fprintf('Statistics (8x8 px):   mean(U) = %+6.4f, mean(V) = %+6.4f, std(U) = %+6.4f, std(V) = %+6.4f\n',...
    meanU2,-meanV2,stdU2,stdV2);
fprintf('Reference:             mean(U) = %+6.4f, mean(V) = %+6.4f, std(U) = %+6.4f, std(V) = %+6.4f\n',...
    8,0,0.7605,0.7751);

%% 
% Compute histogram of u' and show it

% define bin range of histogram
binranges = (5:0.02:11)';

% compute and normalize histogram
bincounts1 = histc(reshape(pivData1.U,N1,1),binranges);
bincounts1 = bincounts1/(sum(bincounts1)*(binranges(2)-binranges(1)));
bincounts2 = histc(reshape(pivData2.U,N2,1),binranges);
bincounts2 = bincounts2/(sum(bincounts2)*(binranges(2)-binranges(1)));

% display it
figure(2);
plot(binranges,bincounts1,'-r',binranges,bincounts2,'-k');
legend('8x8 px 1 pass','8x8 px, 2 passes');
xlabel('displacement U (px)');
ylabel('PDF (a.u.)');
title('Velocity PDF - compare to Fig. 18 in [6]');

%%
% Compute power spectra of u':

% Spectrum from the velocity field after 4 passes
uprime1 = pivData1.U-meanU1;     % velocity difference from the mean velocity
uprime1 = uprime1(1:2:end,:,:);  % reduce amount of velocity data
uspectra1 = zeros(size(uprime1,1)*size(uprime1,3),size(uprime1,2));
% Compute spectrum for every row of velocity data
for ky = 1:size(uprime1,1)
    for kt = 1:size(uprime1,3)
        uspectra1(kt+(ky-1)*size(uprime1,1),:) = abs(fft(uprime1(ky,:,kt))).^2;
    end
end
% compute mean spectrum
uspectra1 = mean(uspectra1,1);
uspectra1 = uspectra1(1:floor(numel(uspectra1)/2));
% Determane wavenumber corresponding to the spectra
dk = 1/(pivData1.X(1,end,1)-pivData1.X(1,1,1));
k1 = (0:(numel(uspectra1)-1))*dk;
% normalize spectrum
uspectra1 = uspectra1/sum(uspectra1) * std(reshape(uprime1,numel(uprime1),1),1)^2/(2*pi*dk) ;

% Spectrum from the velocity field after 5th pass
uprime2 = pivData2.U-meanU2;
uprime2 = uprime2(1:4:end,:,:);  % reduce amount of velocity data
uspectra2 = zeros(size(uprime2,1)*size(uprime2,3),size(uprime2,2));
for ky = 1:size(uprime2,1)
    for kt = 1:size(uprime2,3)
        uspectra2(kt+(ky-1)*size(uprime2,1),:) = abs(fft(uprime2(ky,:,kt))).^2;
    end
end
uspectra2 = mean(uspectra2,1);
uspectra2 = uspectra2(1:floor(numel(uspectra2)/2));
% Determane wavenumber corresponding to the spectra
dk = 1/(pivData2.X(1,end,1)-pivData2.X(1,1,1));
k2 = (0:(numel(uspectra2)-1))*dk;
% normalize spectrum
uspectra2 = uspectra2/sum(uspectra2) * std(reshape(uprime2,numel(uprime2),1),1)^2/(2*pi*dk) ;

% Display the spectra
figure(4);
loglog(2*pi*k1,uspectra1,'-r',2*pi*k2,uspectra2,'-k');
legend('8x8 px 1 pass','8x8 px, 2 passes');
xlabel('k_x (1/px)');
ylabel('E (a.u.)');
xlim([5e-3,1]); ylim([1e-4,100]);
title('Velocity spectra - compare to Figs. 3c and 16 in [6]');

