function [ output ] = DVImgGetMax( Stream,Z,W,T )

if Stream == 0
    helpdlg('Returns the maximum intensity of the section at (Z, W, T).  Note that if operations have been performed on the section since this was last calculated, this value  may be incorrect.','double DVImgGetMax(int Stream,int Z,int W,int T)');
else	
    DVImgCheckZWT(Stream,Z,W,T);
    output = calllib(DVImgLibName,'DVImgGetMax',Stream,Z,W,T);
end
