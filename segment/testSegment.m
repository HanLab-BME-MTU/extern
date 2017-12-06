

% Create an image
image = zeros(256,256);
image(10,:) = 256;
image(50:100,50:100) = 256;
image(150:170,150:170)=256;
%              testSegment.m
%  This code tests the "sbseg" method.  The function "sbseg" must be
%compiled separately before this method is called.
%


image = image+10*randn(256,256);


%  We will use a do-nothing edge detector
edge = ones(256,256);

% Segment with 3 different parameters
u1 = sbseg(image,edge,1e-2);
u2 = sbseg(image,edge,1e-5);
u3 = sbseg(image,edge,2e-6);

% show results
close all;
figure;

subplot(2,2,1);
imagesc(image);
title('original image');

subplot(2,2,2);
imagesc(u1>0.5);
title('mu=1e-2');

subplot(2,2,3);
imagesc(u2>0.5);
title('mu=1e-5');

subplot(2,2,4);
imagesc(u3>0.5);
title('mu=2e-6');
