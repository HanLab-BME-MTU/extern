% Authors:
% Stefan M. Karlsson AND Josef Bigun 
% Please reference the following publication:
% Stefan M. Karlsson AND Josef Bigun, Lip-motion events analysis and lip segmentation using optical flow. CVPR Workshops 2012: 138-145)

% This script calls the function 'vidProcessing'. vidProcessing is our main
% entry point for all video processing. Once its called, a
% figure will open and display the video together with whatever else
% information we desire
% 
% Once you have finished viewing the results of the video processing,
% simply close the window to return back from the function. 

%%%%% argument 'movieType' %%%%%%%%
% This indicates the source of the video to process. You may choose from a
% synthetic video sequence(created on the fly), or load a video through a
% video file (such as 'LipVid.avi'), or to capture video from a connected
% camera(requires the image acquisition toolbox). 
 movieType = 'LipVid.avi'; %assumes a file 'LipVid.avi' in current folder
% movieType = 'camera';
% movieType = 'synthetic';

%%%%% argument 'spdFactor' %%%%%%%%
%%%% a variable to modify the speed of the motion in the synthetic sequence.
%%%% It is ignored unless movieType = 'synthetic'
spdFactor = 1;

%%%%% argument 'method'      %%%%%%%%
%%%%%  optical flow method.  %%%%%%%
% method = 'LK';         %% traditional, explicit Lucas and Kanade
% method = 'LKimproved'; %% regularized parallell Lucas and Kanade (much better and faster than traditional)
% method = 'TS';         %% 3D structure tensor
% method = 'PL';           %% point line flow (the best of the compared)
method = 'PLhsv';        %% point line flow (added components, display by color, high resolution)
% method = 'nothing';    %% output zero fields

%%% There are 2 other options for 'method', both of which                        %%% 
%%% will not give optical flow output, but are useful for testing(and fun):      %%%
% method = 'gradient';  %Displays the gradient values
% method = 'edge';      %Displays the 2D edge detection

%%%%% argument 'bFineScale' %%%%%%%%
%%% determines the scale of differentiation, fine scale otherwise a coarse 
%%% scale. coarse scale gives better stability to large motions, but at the
%%% cost of loosing fine scale information in the video. It determines the
%%% width and height of dx, dy, dt
bFineScale = 1;

%%%%% argument 'nofTimeSlices' %%%%%%%%
%%% the length of temporal integration. This tells how big our sub-volumes
%%% should be for calculations. This gives the depth of dx, dy and dt
%%% This is not to be confused with how many frames are used for creating
%%% a single dx, dy or dt frame. We ALWAYS use only 2 frames for that
nofTimeSlices = 5;

vidProcessing(movieType, method, spdFactor,bFineScale,nofTimeSlices);
