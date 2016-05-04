% Authors:
% Stefan M. Karlsson AND Josef Bigun 

% Please reference the following publication
% Stefan M. Karlsson AND Josef Bigun, Lip-motion events analysis and lip segmentation using optical flow. CVPR Workshops 2012: 138-145)


function make(bDebug)
% This function builds the gradient computation module for you. It will
% attempt to make use of the multi-core parallellism afforded by openMP

if nargin <1
    bDebug = 0;
end

sourceFile = 'Grad3D7do.cpp';

disp('will build the mex-module "Grad3D7do" for your system. Fingers crossed!')

if bDebug
    eval(['mex -g ' sourceFile ' COMPFLAGS="$COMPFLAGS -openmp"  LINKFALGS="$LINKFALGS -openmp"']);
else
    eval(['mex ' sourceFile ' COMPFLAGS="$COMPFLAGS -openmp"  LINKFALGS="$LINKFALGS -openmp"']);
end

disp('Done! Now move the binary "Grad3D7do.mex???"  (where ??? are system specific) to the folder where the original optical flow code is (script runMe.m etc)')
