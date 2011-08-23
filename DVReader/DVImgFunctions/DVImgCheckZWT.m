function [] = DVImgCheckZWT( Stream,Z,W,T )

nz = DVImgGetNumZ(Stream);
nw = DVImgGetNumW(Stream);
nt = DVImgGetNumT(Stream);

if Z < 0 || Z > nz-1
    error('Valid Z range is 0-%d for image open at %d, you called Z of %d.',nz-1,Stream,Z);
elseif W < 0 || W > nw-1
    error('Valid W range is 0-%d for image open at %d, you called W of %d.',nw-1,Stream,W);
elseif T < 0 || T > nt-1
    error('Valid T range is 0-%d for image open at %d, you called T of %d.',nt-1,Stream,T);
end
