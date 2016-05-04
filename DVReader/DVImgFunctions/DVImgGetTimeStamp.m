function [ output ] = DVImgGetTimeStamp( Stream )

if Stream == 0
    helpdlg('Returns the number of seconds between midnight of January 1, 1970 UTC and when the image scan was initiated.  Note that many archived images will not have this data available, as only recent versions of the software began storing it. Images without a time stamp will return 0 when this function is called.','double DVImgGetTimeStamp(int Stream)');
else
    output = calllib(DVImgLibName,'DVImgGetTimeStamp',Stream);
end
