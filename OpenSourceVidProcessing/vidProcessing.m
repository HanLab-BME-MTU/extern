% Authors:
% Stefan M. Karlsson AND Josef Bigun 

% Please reference the following publication:
% Stefan M. Karlsson AND Josef Bigun, Lip-motion events analysis and lip segmentation using optical flow. CVPR Workshops 2012: 138-145)

function [dxOut, dyOut, dtOut, gradInd] = vidProcessing(movieType, method, spdFactor,bFineScale,nofTimeSlices)
% quick usage:
% vidProcessing(); - displays a sequence of test images
% vidProcessing(movieType); - same as above, but on video source indicated by the
% argument. movieType can be:
% 'synthetic' - generates a synthetic video on the fly
% 'camera' - generates video through a connected camera (requires image acquisition toolbox)
% filename - name of video file in the current folder.
% example: vidProcessing('lipVid.avi'); - assumes an avi file is in the current folder
%
% explicit usage:
% [dx,dy,dt] = vidProcessing(movieType);
% Output dx, dy and dt are all WxHxT matrices containing the x, y, and t 
% partial derivatives over time. 
% W and H are the height and width of the video, and T is "nofTimeSlices" 
% used (default T=1). 
% 
% For example, dx(:,:,1), will hold one of the last x derivative images of 
% the video.  
% dx(:,:,1) together with dy(:,:,1) makes up the 2D gradient of one of the
% last images of the video.
%
% [dx,dy,dt] = vidProcessing(movieType, method); displays the same sequence, but
% with a method selected for analyzing the sequence.
% Valid options for 'method' are:
% - "LK"      (Lukas and Kanade)
% - "TS"      (3D structure tensor) 
% - "BONUS"   (User defined algorithm)
% - "NOTHING" ([default] zero fields)
%
% There are 3 special aditional options for 'method', both of which 
% will not give flow output:
% - "gradient"    Displays the gradient values
% - "edge"        Displays the 2D edge detection
% 
% [dx,dy,dt] = vidProcessing('synthetic', method,spdFactor); 
% additionally sets a speed factor (spdFactor) to change the speed of the
% synthetic video generation (if spdFactor=2, then the synthetic sequence 
% is twice as fast)

if nargin <1
    movieType  = 'synthetic'; end
if nargin <2
    method = 'NOTHING'; end
if nargin <3
    spdFactor = 1; end
if nargin <4
    bFineScale = 1; end     
if nargin <5
    nofTimeSlices = 1; end
    
if strcmpi(movieType,'synthetic')
    kindOfMovie  = 'synthetic';
    movieType = '';
    sc =10/spdFactor; %scale vectors for plotting 
elseif strcmpi(movieType,'camera')
    kindOfMovie  = 'camera';
    movieType = '';    
    sc =2; %scale vectors for plotting 
else
	kindOfMovie  = 'file';
    sc =2; %scale vectors for plotting 
end

%a function that sets up the video:
vid =  myVidSetup(kindOfMovie,movieType);

% if running with the camera, now we start it:
if strcmpi(kindOfMovie,'camera')
    start(vid.camIn);
    pause(0.001);
end

%from this point on, we handle the video by the object 'vid'. This is how
%we get the first frame:
curIm = generateFrame(vid, 1,kindOfMovie);

% gradient calculations of the subvolume will be handled by the mex module 
% "Grad3D7do", which is cleared by calling with no args, and initiated by
% calling it with the first frame:
Grad3D7do();
[~,dx, ~, ~, dxC, ~] = Grad3D7do(curIm);

% at this point, "dx" is 2D, and contains the x derivative at fine scale.
% DxC, the coarse scale equivalent. We will work of several frames, so here
% we initiate our datastructures for holding several frames of derivatives:
dx = zeros(size(dx,1),size(dx,2),nofTimeSlices);
dy = dx;
dt = dx;
dxC = zeros(size(dxC,1),size(dxC,2),nofTimeSlices);
dyC = dxC;
dtC = dxC;

if strcmpi(method,'PLhsv')
    edgeIm = zeros(size(curIm,1),size(curIm,2),3);
end

% "bleedOffTerm" is to be used on dx dy dt, so that older values are scaled
% down. This is essentially using a IIR filter in the temporal direction
% The coeffecients of bleedOffTerm are picked so that they will be 100%
% equivalent to a gaussian FIR filtering in the t direction (half a gaussian 
% actually, causal filter). Another good name would be "forgetting factors"
bleedOffTerm = getIIRcoefs(nofTimeSlices);

% index variable for time:
k= 1;

% main loop, runs until user shuts down the figure:
while 1
 k= k+1;
 curIm = generateFrame(vid, k,kindOfMovie,spdFactor);
 % gradInd will be used to index the dx, dy and dt structures (see below)
 gradInd = mod(k,nofTimeSlices)+1;
 
 %implement the bleed off. Gaussian like behaviour:
 for kOffset = 1:(nofTimeSlices-1)
     gradIndOffset = mod(k-kOffset,nofTimeSlices)+1;
     dy(:,:,gradIndOffset)  = dy(:,:,gradIndOffset).*bleedOffTerm(kOffset+1);
     dx(:,:,gradIndOffset)  = dx(:,:,gradIndOffset).*bleedOffTerm(kOffset+1);
     dt(:,:,gradIndOffset)  = dt(:,:,gradIndOffset).*bleedOffTerm(kOffset+1);
	 dyC(:,:,gradIndOffset) = dyC(:,:,gradIndOffset).*bleedOffTerm(kOffset+1);
     dxC(:,:,gradIndOffset) = dxC(:,:,gradIndOffset).*bleedOffTerm(kOffset+1);
     dtC(:,:,gradIndOffset) = dtC(:,:,gradIndOffset).*bleedOffTerm(kOffset+1);
 end
 %calculate the gradient using "Grad3D7do". There are two ways of using the
 %mex module:
%1) in-place computations, using undocumented matlab mex function(faster, but you gotta know what you
% are doing). In versions 2010a and 2011b this has been tested without fault:
%%%%%%%%%%%%%%%%%%% use at own peril:
Grad3D7do(curIm,dy,dx,dt,dyC,dxC,dtC,gradInd);
%%%%%%%%%%%%%%%%%%% 

% %  %stable official Matlab useage(slower):
%  [dy(:,:,gradInd) ,dx(:,:,gradInd), dt(:,:,gradInd),...
%   dyC(:,:,gradInd),dxC(:,:,gradInd),dtC(:,:,gradInd)] = Grad3D7do(curIm);

 if strcmpi(method,'edge')
     if bFineScale
         edgeIm = DoEdgeStrength(dx,dy,gradInd);
     else
         edgeIm = DoEdgeStrength(dxC,dyC,gradInd);
     end
     if k == 2 %if the first run, then setup graphics based on edgeIm
      checkEdgeOutput(edgeIm);
      figH = figure; set(figH, 'Name',method);
      subplot(1,2,1); hImObjEdge = imagesc(edgeIm);
      axis off;axis image;colormap gray(256); title(gca,'Edge Image');
      
      subplot(1,2,2); hImObj = imagesc( curIm,[0,250]); 
      axis off;axis image;colormap gray(256);title(gca,'original sequence');
     end
 elseif strcmpi(method,'PLhsv')
     if bFineScale
    	[U1, V1] = DoFlowPL(dx,dy,dt);
     else
        [U1, V1] = DoFlowPL(dxC,dyC,dtC);
     end
     edgeIm(:,:,1) = (atan2(imresizeNN(V1,size(curIm)),imresizeNN(U1,size(curIm)))+ pi)/(2*pi);          
     edgeIm(:,:,2) = min(1,imresizeNN(U1.^2 + V1.^2,size(curIm)));
     edgeIm(:,:,3) = double(curIm)/255;
     if k == 2 %if the first run, then setup graphics based on edgeIm
      figH = figure; set(figH, 'Name',method);
      hImObjEdge = image(hsv2rgb(edgeIm));
      axis off;axis image; title(gca,'PLhsv Image');
    end             
 elseif strcmpi(method,'gradient')
    if k == 2 %if the first run, then setup graphics based on dx, dy and dt
      figH = figure; set(figH,'Name','Displaying 3D gradient');
      plotRange = 0.08;
      if bFineScale
          subplot(2,2,1); hImObjDx = imagesc(dx(:,:,gradInd),[-plotRange,plotRange]); 
          axis off;axis image;colormap gray(256);title(gca,'dx');

          subplot(2,2,2); hImObjDy = imagesc(dy(:,:,gradInd),[-plotRange,plotRange]); 
          axis off;axis image;colormap gray(256);title(gca,'dy');

          subplot(2,2,3); hImObjDt = imagesc(dt(:,:,gradInd),[-plotRange,plotRange]); 
          axis off;axis image;colormap gray(256);title(gca,'dt');
      else
          subplot(2,2,1); hImObjDx = imagesc(dxC(:,:,gradInd),[-plotRange,plotRange]); 
          axis off;axis image;colormap gray(256);title(gca,'dxC');

          subplot(2,2,2); hImObjDy = imagesc(dyC(:,:,gradInd),[-plotRange,plotRange]); 
          axis off;axis image;colormap gray(256);title(gca,'dyC');

          subplot(2,2,3); hImObjDt = imagesc(dtC(:,:,gradInd),[-plotRange,plotRange]); 
          axis off;axis image;colormap gray(256);title(gca,'dtC');
      end      
	  subplot(2,2,4); hImObj = imagesc( curIm,[0,250]); 
      axis off;axis image;colormap gray(256);title(gca,'original sequence');
    end
 else %flow methods
   if bFineScale
    [U1, V1, U2, V2] = DoFlow(dx,dy,dt,method);
   else
    [U1, V1, U2, V2] = DoFlow(dxC,dyC,dtC,method);
   end
    if k == 2 %if the first run, then setup graphics based on size of flow field
      checkFlowOutput(U1, V1, U2, V2);
      flowRes = max(size(U1));
      [hImObj,hQvObjPoints,hQvObjLines] =  myGraphicsSetup(vid,kindOfMovie,flowRes);
      title(gca,['current image. Flow vectors are scaled by ' num2str(sc) ':1' ]);
      figH = gcf; set(figH, 'Name',method);
    end
 end
 
 % UPDATE GRAPHICS
 if ishandle(figH)%check if the figure is still open
    if strcmpi(method,'edge')
        set(hImObj ,'cdata',curIm);
        set(hImObjEdge ,'cdata',edgeIm);
    elseif strcmpi(method,'PLhsv')
        set(hImObjEdge ,'cdata',hsv2rgb(edgeIm));        
    elseif strcmpi(method,'gradient')
        if bFineScale
            set(hImObjDx,'cdata',dx(:,:,gradInd));
            set(hImObjDy,'cdata',dy(:,:,gradInd));
            set(hImObjDt,'cdata',dt(:,:,gradInd));
        else
            set(hImObjDx,'cdata',dxC(:,:,gradInd));
            set(hImObjDy,'cdata',dyC(:,:,gradInd));
            set(hImObjDt,'cdata',dtC(:,:,gradInd));
        end
            set(hImObj  ,'cdata',curIm);
    else %flow methods
        set(hImObj ,'cdata',curIm);
        set(hQvObjLines ,'UData', sc*U1, 'VData', sc*V1);
        set(hQvObjPoints,'UData', sc*U2, 'VData', sc*V2);
    end
 else%user has killed the main figure, break the loop and exit
    break;
 end

% a brief pause lets matlab update figures, and is good for computer 
% stability. Change this value to introduce lag-time.
pause(0.001);
end

% Clears up memory:
Grad3D7do();
if strcmpi(kindOfMovie,'file')
    delete(vid);   
elseif strcmpi(kindOfMovie,'camera')
	delete(vid.camIn);   
end

if bFineScale
    dxOut = dx;dyOut = dy; dtOut = dt;
else
    dxOut = dxC;dyOut = dyC; dtOut = dtC;
end



%%%%% THE REST OF THIS FILE CONTAINS ONLY HELPER FUNCTIONS
% myVidSetup- setups the video feed. 3 types are handled different,
% depending on the 'kindOfMovie' input arg. Can be 'file'(file on disk),
% 'synthetic'(manufactured test sequence) or 'camera' (setups the default 
% video input device for capturing video for this application)
function vid =  myVidSetup(kindOfMovie,movieType)
if strcmpi(kindOfMovie, 'file')
        %we open the file for reading
        vid = mmreader(movieType);
    elseif strcmpi(kindOfMovie, 'synthetic') 
        %this option just requires two fields in the 'vid' structure
        vid.Height = 128;   vid.Width  = 128;
    elseif strcmpi(kindOfMovie, 'camera')
        vid.Height = 128;   vid.Width  = 128;
%         first, reset image aqcuisition devices. This also tests if the
%         toolbox is available. If not, exit with error
	try      
        imaqreset;
    catch %#ok<CTCH>
        error('Image Aquisition toolbox not available. Camera unsupported');
	end

%     get info on supported formats for this capture device
    dev_info = imaqhwinfo('winvideo',1);
    strVid = dev_info.SupportedFormats;
    
%     strVid contains all the supported formats as strings , for example 
%     'I420_320x240' is one such format. I want to pick the smallest
%     resolution available, so I parse these strings in what follows
    splitStr = regexpi(strVid,'x|_','split');
    pickedFormat = 0;
    resolutionFormat = Inf;
    for ik = 1:length(strVid)
        resW = str2double(splitStr{ik}{2});
        resH = str2double(splitStr{ik}{3});
        if (resW > (vid.Width-1) )&&(resH > (vid.Height-1) )&& (resW*resH)<resolutionFormat
            resolutionFormat = (resW*resH);
            pickedFormat = ik;
        end
    end
    % pick the selected format, color and a region of interest 128x128 big:
    vid.camIn = videoinput('winvideo',1,strVid{pickedFormat});
    set(vid.camIn, 'ReturnedColorSpace', 'gray');
    set(vid.camIn, 'ROIPosition', [1 1 vid.Width vid.Height]);
    %let the video go on forever, grab one frame every update time:
    triggerconfig(vid.camIn, 'manual');
    src = getselectedsource(vid.camIn);
    frameRates = set(src, 'FrameRate');
    src.FrameRate = frameRates{1};
    
end
    
function [hImObj,hQvObjPoints,hQvObjLines] =  myGraphicsSetup(vid,kindOfMovie,flowRes)

    imSize = vid.Height;
    figure;
    curIm = generateFrame(vid, 1,kindOfMovie);
    hImObj   = imagesc( curIm,[0,250]);
    colormap gray;axis off;axis manual;axis image;
    hold on;
    axisInterval = linspace(1,imSize,flowRes+2);
    axisInterval = axisInterval(2:end-1) ;

    hQvObjLines = quiver(axisInterval,axisInterval, zeros(flowRes),  zeros(flowRes),0 ,'m','MaxHeadSize',5,'Color',[.9 .2 .1]);%, 'LineWidth', 1);
    hQvObjPoints= quiver(axisInterval,axisInterval, zeros(flowRes),  zeros(flowRes),0 ,'m','MaxHeadSize',5,'Color',[.1 .9 .2]);%, 'LineWidth', 1);    
    truesize([imSize imSize]*4);
    
function newIm = generateFrame(vid, t,kindOfMovie,spdFactor)
persistent iix iiy;
if nargin<4
    spdFactor = 1;end

if strcmp(kindOfMovie,'file') 
    %get frame at t in video file, where t is "bouncing" forward and
    %backward (thus, video is played for infinity, by playing itself backwards half of the time)
    t=varyT(round(t),vid.NumberOfFrames);
    newIm = rgb2gray(read(vid, t));
elseif strcmpi(kindOfMovie,'synthetic')  % generate test image:
    %%%%%%%%%%% image generation Params %%%%%%%%%%%
    spd    = spdFactor/300; %speed of motion of the patterns generated
    cen1   = 0.45;          %centre of circle 1
    cen2   = -cen1;         %centre of circle 2
    cen3   = -cen1;         %centre of circle 3, offset in y
    cW     = 0.27;          %radius of circles
    cFuz   = 2.7;           %fuzziness of the boundary
    cDetail= 0.65;          %detail of the interior pattern(frequency of sinusoids)

    if isempty(iix)
        [iix,iiy] = meshgrid(linspace(-1,1,vid.Height));
    end
    iX = (iix    -0.2*sin(pi*(t+0.5)*spd*3*2));
    iY = (iiy-0.22-0.45*cos(pi*(t+0.5)*spd*2  ));

   newIm = uint8(... 
     185*(0.4+cos(2+iY*cDetail*16*pi).*cos(2+iX*cDetail*16*pi)).*sig(cW^2-(iX+cen1).^2-iY.^2, cFuz*cW/50) + ... //disk 1
     128*(1+                             cos(iX*cDetail*20*pi)).*sig(cW^2-(iX+cen2).^2-iY.^2, cFuz*cW/50) + ... //disk 2
     128*(1+cos(iY*cDetail*20*pi)                             ).*sig(cW^2-(iY-cen3).^2-iX.^2, cFuz*cW/50));  ... //disk 3
else %if strcmpi(kindOfMovie,'camera')  % capture from camera:
    newIm = fliplr(squeeze(getsnapshot(vid.camIn)));
end

% checks for consistency in the DoFlow output
function checkFlowOutput(U1, V1, U2, V2)
      if (ndims(U1) ~=2 || ndims(U2) ~=2|| ndims(V1) ~=2|| ndims(V2) ~=2)
          error('function "DoFlow" has returned invalid output. "U1", "V1", "U2" and "V2" must all be 2D');end      
      if (size(U1,1) ~= size(U1,2))
          error('function "DoFlow" has returned an invalid output. "U1", "V1", "U2" and "V2" must all be square (height=width)');end
      if (sum(size(U1) ~= size(U2))||sum(size(U1) ~= size(V1))||sum(size(U1) ~= size(V2)))
          error('function "DoFlow" has returned an invalid output. "U1", "V1", "U2" and "V2" must all be of equal size');end
      
% checks for consistency in the DoEdgeStrength output
function checkEdgeOutput(edgeIm)
      if (ndims(edgeIm) ~=2 || size(edgeIm,1) ~= size(edgeIm,2))
          error('function "DoEdgeStrength" has returned an invalid edge image. edgeIm should be square, 2D');
      end      

%%%%%%% a sigmoidal function. It does the logical operation (x>0) in a fuzzy manner %%%%
function out=sig(x,fuzziness)
if nargin <2
    fuzziness = 1;end;
if (fuzziness == 0)
    out = x > 0;
else
    out= (1+erf(x./fuzziness))/2;
end

%generate coefficients for IIR filtering in temporal direction of the
%gradient image
function IIRcoeffs = getIIRcoefs(n)

g= gaussgen(0.5*(n-1),(n-1)*2+1);
g=g(n:end)/g(n);
IIRcoeffs = ([1 g(2:n)./g(1:n-1)]).^(1/2);

%%% varies t to "bounce-loop" a video of length T. It goes on forever, by
%%% alternating playing the video forwards and backwards.
function t2 = varyT(t,T)
    T = T-1;
    t= t-1;
    t2 = mod(t,T)-2*mod(floor((t)/T),2).*(mod(t,T)-T/2)+1;
    
function h=gaussgen(std,siz)
% This function generates a 1-D gaussian kernel
x2 = (-(siz-1)/2:(siz-1)/2).^2;
h = exp(-(x2)/(2*std*std));
h = h/sum(sum(h));

function outputImage = imresizeNN(inputImage, newSize)
%%%%%%% imresizeNN(inputImage, newSize) is identical to built in 
%%%%%%% imresize(inputImage, newSize, 'nearest'), but is much faster
oldSize = size(inputImage);  
scale = newSize./oldSize;    

% Compute a resampled set of indices:
rowIndex = min(round(((1:newSize(1))-0.5)./scale(2)+0.5),oldSize(1));
colIndex = min(round(((1:newSize(2))-0.5)./scale(2)+0.5),oldSize(2));

outputImage = inputImage(rowIndex,colIndex);
