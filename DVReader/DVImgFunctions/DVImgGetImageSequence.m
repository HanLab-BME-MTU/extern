function [ output ] = DVImgGetImageSequence( Stream )

if Stream == 0
    helpdlg('Returns the order that the XY sections are stored in on the storage device.  This can be 0 (ZTW), 1 (WZT), or 2 (ZWT).  2, ZWT, is a common sequence.  Although a well-written script can handle any image sequence, it can be advantageous, time-wise, to process the sections in the same order they are stored to keep from skipping around on the storage device.','int DVImgGetImageSequence(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetImageSequence',Stream);
end
