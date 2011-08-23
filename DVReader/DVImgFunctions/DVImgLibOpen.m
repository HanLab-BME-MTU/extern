function [] = DVImgLibOpen(ListFunctions)
if nargin < 1
    ListFunctions = 0;
end

LibPath = dvimgroot();
tf=ispc;
if tf == 1
   SharedPath = fullfile(LibPath, [DVImgLibName '.dll']);
else
   SharedPath = fullfile(LibPath, [DVImgLibName '.so']); 
end
HeaderPath = fullfile(LibPath, [DVImgLibName '.h']);
MPath = fullfile(LibPath, [DVImgLibName]);

ThisDir = pwd;
cd(LibPath);
loadlibrary(SharedPath,@dv_img_io_m);
cd(ThisDir);

if ListFunctions ~= 0
    libfunctions (DVImgLibName) -full
end
