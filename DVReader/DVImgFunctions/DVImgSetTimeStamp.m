function [ output ] = DVImgSetTimeStamp( Stream,TimeStamp );

if Stream == 0
    helpdlg('Sets the time stamp of the requested stream.  TimeStampSecs should be in the format of the number of seconds since midnight of Jan 1, 1970, UTC.  The first valid time-stamp is 315561600 secs, which corresponds to Jan 1, 1980.  Not all DV images can accept a time stamp.','int DVImgSetTimeStamp(int Stream,double TimeStampSecs)');
else
    output = calllib(DVImgLibName,'DVImgSetTimeStamp',Stream,TimeStamp);
    DVImgPrintErrText(output);
end
