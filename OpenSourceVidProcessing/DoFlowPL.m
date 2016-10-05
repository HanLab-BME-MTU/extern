% Authors:
% Stefan M. Karlsson AND Josef Bigun 

% Please reference the following publication
% Stefan M. Karlsson AND Josef Bigun, Lip-motion events analysis and lip segmentation using optical flow. CVPR Workshops 2012: 138-145)

function [U1, V1, U2, V2] = DoFlowPL(dx,dy,dt)
% function DoFlow inputs 3D volume images, dx, dy, dt, corresponding to the
% 3D gradients in a spatio-temporal image volume.
% the dimensions are
% [height, width, T] =size(dx)=size(dy)=size(dt)
% where T is the time interval in which the video is considered for optical
% flow calculation
% 
% The output, [U1, V1, U2, V2] are the component images of point motion
% (U1, V1) and line motion (U2, V2). If the function is called with only
% two outputs as:
%     [U, V] = DoFlowPL(dx,dy,dt)
% then U = U1+U2, V = V1+V2. This will be quite similar to a locally
% regularized Lucas Kanade, but better results and no need for a
% regularization parameter

%     MOMENT CALCULATIONS
gaussStd = 1.5; % for tensor smoothing.
gg  = gaussgen(gaussStd); %% filter for tensor smoothing

m200=filter2(gg',filter2(gg,sum(dx.^2 ,3)))/size(dx,3);
m020=filter2(gg',filter2(gg,sum(dy.^2 ,3)))/size(dx,3);
m110=filter2(gg',filter2(gg,sum(dx.*dy,3)))/size(dx,3);
m101=filter2(gg',filter2(gg,sum(dx.*dt,3)))/size(dx,3);
m011=filter2(gg',filter2(gg,sum(dy.*dt,3)))/size(dx,3);

% not really a regularization paramater, just to avoid division by zero:
TikConst  = 0.00000001; 

tensorTr = (m200 + m020+TikConst);
if nargout == 4
    V1 = -(((m200-m020).^2 + 4*m110.^2).*m011)./tensorTr.^3;
    U1 = -(((m200-m020).^2 + 4*m110.^2).*m101)./tensorTr.^3;
    V2 = 4*( m101.*m110 - m011.*m200)./tensorTr.^2;
    U2 = 4*(-m101.*m020 + m011.*m110)./tensorTr.^2;
else
% If one does not want the seperation into point and lines, one can simply add the
% components together and get superior local regularization:
    V1 = 4*( m101.*m110 - m011.*m200)./tensorTr.^2 ...
         -(((m200-m020).^2 + 4*m110.^2).*m011)./tensorTr.^3;
    U1 = 4*(-m101.*m020 + m011.*m110)./tensorTr.^2 ...
        -(((m200-m020).^2 + 4*m110.^2).*m101)./tensorTr.^3;
end

%%%% the rest of this file contains helper functions
function h=gaussgen(std)
% This function generates a 1-D gaussian kernel
siz = round(5*std);
siz =siz + mod(siz-1,2); %make sure odd sized

x2 = (-(siz-1)/2:(siz-1)/2).^2;
h = exp(-(x2)/(2*std*std));
h = h/sum(sum(h));
