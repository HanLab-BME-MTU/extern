function mask = circleMask(mskSize)
%CIRCLEMASK creates a nd-circular or ellipsoidal mask for use with e.g. morphological operators 
%
% SYNOPSIS: mskSize = mask(mskSize)
%
% INPUT 		mskSize : 1-by-nDims array with the size of the mask in
%                         each dimension.
%
% OUTPUT        mask : logical array of size mskSize with ones for each
%                      voxel that is inside a 'n-d circle' with radius
%                      floor(mskSize/2)
%
% REMARKS       
%
% SEE ALSO strel
%
% created with MATLAB ver.: 7.10.0.499 (R2010a) on Microsoft Windows 7 Version 6.1 (Build 7600)
%
% created by: Jacques Boisvert
% DATE: 04-Jun-2010
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Validation

if nargin < 1 || isempty(mskSize)
    error('You need to provide the mask size so that a mask can be created');
end

%% Circle Mask

%Number of dimension
nDim = length(mskSize);
range = cell(nDim,1);

%I calculate the range of my coordinates. Let them go from -1 to 1 so that
%thresholding becomes a lot easier below. Basically, we're doing a
%coordinate transformation here to make it look like we're creating a
%circle, even though the mask may be ellipsoidal.
for idx = 1:nDim
    range{idx} = linspace(-1,1,mskSize(idx));
end
% Create a coordinate grid for each dimension
[coordinates{1:nDim}] = ndgrid(range{:});
tmp = cat(nDim+1,coordinates{:});
%I keep the one that are inside the radius of 1.
mask = sum(tmp.^2 , nDim+1) <= 1;