% Authors:
% Stefan M. Karlsson AND Josef Bigun 

% Please reference the following publication
% Stefan M. Karlsson AND Josef Bigun, Lip-motion events analysis and lip segmentation using optical flow. CVPR Workshops 2012: 138-145)

function [U1, V1, U2, V2] = DoFlow(dx,dy,dt,method)
% function DoFlow inputs 3D volume images, dx, dy, dt, corresponding to the
% 3D gradients in a spatio-temporal image volume.
% the dimensions are
% [height, width, T] =size(dx)=size(dy)=size(dt)
% where T is the time interval in which the video is considered for optical
% flow calculation
% 
% The output, [U1, V1, U2, V2] are the component images of point motion
% (U1, V1) and line motion (U2, V2)

if nargin<4
%%% possibilities for 'method' are:  %%
%     method = 'LK';         %% Lucas Kanade
%     method = 'LKimproved'; %% Lucas Kanade improved
%     method = 'TS';         %% 3D structure tensor
%     method = 'PL';         %% Point-line flow (Stefan M. Karlsson AND Josef Bigun, Lip-motion events analysis and lip segmentation using optical flow. CVPR Workshops 2012: 138-145)
    method = 'NOTHING';      %% Do nothing, i.e. return zeros
end

flowRes = 21; %resolution of the flow field in the image

%     MOMENT CALCULATIONS
      gaussStd = 1.5; % for tensor smoothing.
      gg  = gaussgen(gaussStd); %% filter for tensor smoothing

%     moment m200, calculated in three steps explicitly
%     1) make elementwise product, and sum along time direction (time integration):
      momentIm = sum(dx.^2,3);
     
%     2) smooth with large seperable gaussian filter (spatial integration)
      momentIm = filter2(gg',filter2(gg,momentIm));

%     3) downsample to specified resolution:     
      m200 =  imresizeNN(momentIm ,[flowRes flowRes])/size(dx,3);

%    The remaining moments are calculated in EXACTLY the same way as above, condensed to one liners:
 m020=imresizeNN(filter2(gg',filter2(gg,(sum(dy.^2 ,3)))),[flowRes flowRes])/size(dx,3);
 m002=imresizeNN(filter2(gg',filter2(gg,(sum(dt.^2 ,3)))),[flowRes flowRes])/size(dx,3);
 m110=imresizeNN(filter2(gg',filter2(gg,(sum(dx.*dy,3)))),[flowRes flowRes])/size(dx,3);
 m101=imresizeNN(filter2(gg',filter2(gg,(sum(dx.*dt,3)))),[flowRes flowRes])/size(dx,3);
 m011=imresizeNN(filter2(gg',filter2(gg,(sum(dy.*dt,3)))),[flowRes flowRes])/size(dx,3);

 % Thresholds, and regularization paramaters, play around freely:

EPSILONLK = 0.0000001;  
EPSILONTS = 0.00001; GAMMA=0.333;
TikConst  = 0.00001; 

if strcmpi(method,'LK')
%initialize the two fields to zero:
 U1 = zeros(size(m200));
 V1 = U1;

%%%%%%%%%%%%%% flow by the Lucas and Kanade algorithm  %%%%%%%%%% 
%%%%% this is a traditional, explicit way of estimating flow. It suffers
%%%%% from the problems of numerical instability and slow implementation.
%%%%% never use this implementation for anything but educational purposes
%%%%% instead, see the next method "LKimproved"
  for x=1:size(m011,1)
    for y=1:size(m011,2)
        %%%build the 2D structure tensor
        S2D  =... 
            [m200(x,y), m110(x,y);...
             m110(x,y), m020(x,y)];

        %%%%check that S2D is invertible :
        if(det(S2D)>EPSILONLK)
            b = [m101(x,y);...
                 m011(x,y)];
            %%%% calculate the velocity vector by the relation 
            %%%% between vector b, and matrix S2D (2D structure tensor)
            v = -S2D\b;
            U1(x,y) = v(1);
            V1(x,y) = v(2);
        end
    end
  end
elseif strcmpi(method,'LKimproved')
%%%%%%%%%%%%%%% Parallell and regularized Lukas Kanade: %%%%%%%%%%     
%%%% this introduces stability to linear motions. It does not let
%%%% you label the motions as 'point' or 'line'. This method is superior to
%%%% the traditional implementation of Lucas annd Kanade algorithm in every
%%%% way (previous method is kept for educational purposes). Even though we
%%%% do not explicitly deal with a matrix inverse, it is implicit. You can
%%%% easily derive the expressions below using symbolic programming.
%%%% E.g, with Matlab symbolic toolbox you can write:
%          syms m200 m020 m110 m101 m011;
%          b   = [m101; m011];
%          S2D = [m200, m110; ...
%                 m110, m020];
%          v = -S2D\b;
%%%% this would derive the expressions below for the estimated motion

m200 = (m200 + TikConst);
m020 = (m020 + TikConst);
tensDet = (m020.*m200 - m110.^2);
U1 = (-m101.*m020 + m011.*m110)./tensDet;
V1 =( m101.*m110 - m011.*m200)./tensDet;
  
elseif strcmpi(method,'TS')
%%%%%%%%%%%%%% flow by the 3D structure tensor algorithm  %%%%%%%%%% 
%%%%% This is the traditional 3D structure tensor algorithm, that suffers
%%%%% from some similar issues as the traditional LK algorithm. It can,
%%%%% however, distinguis automatically between line motions and point
%%%%% motions. The method 'PL' below is an improvement to this in most
%%%%% cases.
 %initialize the two fields to zero:
 U1 = zeros(size(m200));
 V1 = U1;
 U2 = U1;
 V2 = U1;

for x=1:size(m011,1)
    for y=1:size(m011,2)
    %%%build the 3D structure tensor
    S3D =  [m200(x,y), m110(x,y), m101(x,y);...
            m110(x,y), m020(x,y), m011(x,y);...
            m101(x,y), m011(x,y), m002(x,y)];

    % Compute the eigenvalues
	[eigenvects, eigenvals]=eig(S3D);
	eigenvals=diag(eigenvals);

	% And sort them. eigenvals(1) is the smallest; eigenvals(3)
	% is the largest. 
	[eigenvals, index]=sort(eigenvals);   

   % eigenvals=log(eigenvals);
    sumeig=sum(eigenvals);
    %%% We only proceed if we are in a spatio temporal region with
    %%% variation. The sum of the eigenvalues tell us this
    %%% (if we are in a constant region, sumeig=0 and we skip it):
    if abs(sumeig) >EPSILONTS
        CER_L=(eigenvals(3)-eigenvals(2))/eigenvals(3);
        CER_P=(eigenvals(2)-eigenvals(1))/eigenvals(3);
        CER_B=(eigenvals(1)             )/eigenvals(3);
        [CER IND]=max([CER_L CER_P CER_B]);

        if  CER > GAMMA
            switch IND
                case 1 %Line motion
                    % In this case we have spatial orientation and constant motion.
                    % Only one eigenvalue is different from zero;
                    % the corresponding eigenvector gives us the normal velocity
                    e=eigenvects(:,index(3)); % largest eigenvalue
                    denom=(e(1)^2+e(2)^2);
                    if denom>EPSILONTS
                        U1(x,y)=- e(3)*e(1)/denom;
                        V1(x,y)=- e(3)*e(2)/denom;
                    end
                case 2 %Point motion
                    % Case of distributed spatial structure and constant motion
                    e=eigenvects(:,index(1)); % smallest eigenvalue
                    if e(3)>EPSILONTS
                        U2(x,y)=e(1)/(e(3));
                        V2(x,y)=e(2)/(e(3));
                    end
                otherwise
                    % Distributed spatial structure and non-constant motion. We do nothing in this case
            end
        end
     end %if sumeig
    end % for x
end
 
 elseif strcmpi(method,'PL') % point line flow, from:
    % Stefan M. Karlsson AND Josef Bigun, Lip-motion events analysis and lip segmentation using optical flow. CVPR Workshops 2012: 138-145)
    % This seems the best method of those implemented here. 

    tensorTr = (m200 + m020+TikConst);

    if nargout == 4
        V1 = -(((m200-m020).^2 + 4*m110.^2).*m011)./tensorTr.^3;
        U1 = -(((m200-m020).^2 + 4*m110.^2).*m101)./tensorTr.^3;
        V2 = 4*( m101.*m110 - m011.*m200)./tensorTr.^2;
        U2 = 4*(-m101.*m020 + m011.*m110)./tensorTr.^2;
    else
    % If one does not want the seperation into point and lines, one can simply add the
    % components together and get superior regularization to any of the
    % above methods by:
        V1 = 4*( m101.*m110 - m011.*m200)./tensorTr.^2 ...
             -(((m200-m020).^2 + 4*m110.^2).*m011)./tensorTr.^3;
        U1 = 4*(-m101.*m020 + m011.*m110)./tensorTr.^2 ...
            -(((m200-m020).^2 + 4*m110.^2).*m101)./tensorTr.^3;
    end

elseif strcmpi(method,'NOTHING')
    return;
else %then a method string we dont recognize => Error
    error(['undefined method "' method '". Choose from "LK", "LKimproved", "TS", "PL", or "NOTHING".']);
end


%%%% the rest of this file contains helper functions

function outputImage = imresizeNN(inputImage, newSize)
%%%%%%% imresizeNN(inputImage, newSize) is identical to built in 
%%%%%%% imresize(inputImage, newSize, 'nearest'), but is much faster
oldSize = size(inputImage);  
scale = newSize./oldSize;    

% Compute a resampled set of indices:
rowIndex = min(round(((1:newSize(1))-0.5)./scale(2)+0.5),oldSize(1));
colIndex = min(round(((1:newSize(2))-0.5)./scale(2)+0.5),oldSize(2));

outputImage = inputImage(rowIndex,colIndex);


function h=gaussgen(std)
% This function generates a 1-D gaussian kernel
siz = round(5*std);
siz =siz + mod(siz-1,2); %make sure odd sized

x2 = (-(siz-1)/2:(siz-1)/2).^2;
h = exp(-(x2)/(2*std*std));
h = h/sum(sum(h));
