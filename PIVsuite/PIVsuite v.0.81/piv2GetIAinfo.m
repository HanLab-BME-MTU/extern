function [IAinfo] = piv2GetIAinfo(imSize,pivPar)
% pivGetIAinfo - create an array containing information about interrogation areas
%
% Inputs:
%    imSize ... vector of two components with number of rows and columns in image
%    pivPar ... structure containing these fields:
%       iaSizeX, iaSizeY ... size of interrogation areas 
%       iaStepX, iaStepY ... grid spacing between IA's
%
% Outputs:
%    IAinfo ... 2D array of structures, containing information about IA. Structure has these fiels:
%       X, Y ... centers of IA's
%       minX,maxX ... minimum and maximum X coordinate of pixels belonging to the IA
%       minY,maxY ... minimum and maximum Y coordinate of pixels belonging to the IA

%% 0. Extract from pivPar used fields (for shortenning the code);
iaSizeX = pivPar.iaSizeX;
iaSizeY = pivPar.iaSizeY;
iaStepX = pivPar.iaStepX;
iaStepY = pivPar.iaStepY;
imSizeX = imSize(2);
imSizeY = imSize(1);

%% 1. Get the number of IA's
iaNX = floor((imSizeX - iaSizeX)/iaStepX)+1;
iaNY = floor((imSizeY - iaSizeY)/iaStepY)+1;

%% 2. Distribute IA's (undeformed image, no offset):
auxLengthX = iaStepX * (iaNX-1) + iaSizeX;
auxLengthY = iaStepY * (iaNY-1) + iaSizeY;
auxFirstIAX = floor((imSizeX - auxLengthX)/2) + 1;
auxFirstIAY = floor((imSizeY - auxLengthY)/2) + 1;
iaStartX = (auxFirstIAX:iaStepX:(auxFirstIAX+(iaNX-1)*iaStepX))';  % first columns of IA's
iaStartY = (auxFirstIAY:iaStepY:(auxFirstIAY+(iaNY-1)*iaStepY))';  % first rows of IA's
iaStopX = iaStartX + iaSizeX - 1;  % last columns of IA's
iaStopY = iaStartY + iaSizeY - 1;  % last rows of IA's
iaCenterX = (iaStartX + iaStopX)/2;    % center of IA's (usually between pixels)
iaCenterY = (iaStartY + iaStopY)/2;

%% 3. Create output structure
[X,Y] = meshgrid(iaCenterX,iaCenterY);
[minX,minY] = meshgrid(iaStartX,iaStartY);
[maxX,maxY] = meshgrid(iaStopX,iaStopY);

for ky=1:size(X,1)
    for kx=1:size(X,2)
        IAinfo(ky,kx).X = X(ky,kx);          %#ok<*AGROW>
        IAinfo(ky,kx).Y = Y(ky,kx);
        IAinfo(ky,kx).minX = minX(ky,kx);
        IAinfo(ky,kx).minY = minY(ky,kx);
        IAinfo(ky,kx).maxX = maxX(ky,kx);
        IAinfo(ky,kx).maxY = maxY(ky,kx);
    end
end
