/* ==================================================================*/
/*
   dv_img_io.h

   Public include file for the DV Image IO Library

   Version    Date       Author  Comments
   1.00       02/25/07   CSB     First working version ready for testing.

---------------------------------------------------------------------
(C) Copyright - 2007 - Applied Precision LLC
CONFIDENTIAL - This document and the information contained herein are the
confidential and proprietary property of Applied Precision, LLC, and may
not be reproduced, disclosed to any person, or used for any purpose
whatsoever without the prior written approval of Applied Precision, LLC
*/

/* ==================================================================*/

#ifdef __cplusplus
extern "C" {
#endif

#ifdef WIN32
  // needed for external .def file used to export functions in the DLL
  #undef AFX_DATA
  #define AFX_DATA AFX_EXT_DATA
#endif

// Basic Information and Variable Conventions:
// -------------------------------------------
// Pixel sizes and XYZ coordinates are in microns.
// Wavelengths are in nanometers.
// Wavelength, z-section, and time-point numbers are 0 based.
// Image stream numbers should be between 1 and 20.
// The display scaling exponent "dScaleExp" used in DVImgSetDisplay() is typically 1.
// Stream attributes are either "new" or "old"
//
// Image sequence code numbers. [See "iImgSequence" in the prototype for DVImgCreate()]
//   0    Z -> T -> W
//   1    W -> Z -> T
//   2    Z -> W -> T
//
// Image pixel data type numbers.  [See "iPixDataType" in the prototype for DVImgCreate()]
//   0    8-bit integer
//   1    16-bit, signed integer
//   2    32-bit floating-point
//   3    Unused/unsupported.
//   4    Unused/unsupported.
//   5    Unused/unsupported.
//   6    16-bit, unsigned integer
//   7    32-bit, signed integer
//   8    64-bit floating-point

//---------------------------------------------------------------------
// Public Function Prototypes                                                 

// Basic image access functions.
int    DVImgOpen(int iStream, char *cFileName, char *cAccessType);
int    DVImgCreate(int iStream, char *cFileName,
                         int iNX,int iNY,int iNZ,int iNW,int iNT,
                         double dPixSizeX,double dPixSizeY,double dPixSizeZ,
                         int iImgSequence,int iPixDataType,int iExtHeaderType);
int    DVImgClose(int iStream);

// Image parameter alteration functions
int    DVImgSetDisplay(int iStream,int iWaveNum,double dScaleMin, double dScaleMax, double dScaleExp);
int    DVImgSetIntenStats(int iStream, int iWaveNum,double dMin,double dMax,double dMean);
int    DVImgSetOrigin(int iStream,double dXOrig,double dYOrig,double dZOrig);
int    DVImgSetPixelSize(int iStream, double dDX, double dDY, double dDZ);
int    DVImgSetPosX(int iStream,int iZ,int iW,int iT,double dPosX);
int    DVImgSetPosY(int iStream,int iZ,int iW,int iT,double dPosY);
int    DVImgSetPosZ(int iStream,int iZ,int iW,int iT,double dPosZ);
int    DVImgSetTime(int iStream,int iZ,int iW,int iT,double dT);
int    DVImgSetPhotoSen(int iStream,int iZ,int iW,int iT,double dPhotoVal);
int    DVImgSetMin(int iStream,int iZ,int iW,int iT,double dMin);
int    DVImgSetMax(int iStream,int iZ,int iW,int iT,double dMax);
int    DVImgSetMean(int iStream,int iZ,int iW,int iT,double dMean);
int    DVImgSetTimeStamp(int iStream,double dTimeStampSecs);
int    DVImgSetVerboseMode(int iVerboseModeOn);
int    DVImgSetLensID(int iStream,int iLensID);
int    DVImgSetWavelength(int iStream, int iWaveNum, double dWavelength);
int    DVImgSetTypeConversionState(int iStream,int iDataTypeConversionFlag);

// Query functions
int    DVImgGetBytesPerPixel(int iStream);
int    DVImgGetDataType(int iStream);
int    DVImgGetImageSequence(int iStream);
double DVImgGetDataTypeMin(int iStream);
double DVImgGetDataTypeMax(int iStream);
double DVImgGetPixelSizeX(int iStream);
double DVImgGetPixelSizeY(int iStream);
double DVImgGetPixelSizeZ(int iStream);

double DVImgGetPosX(int iStream,int iZ,int iW,int iT);
double DVImgGetPosY(int iStream,int iZ,int iW,int iT);
double DVImgGetPosZ(int iStream,int iZ,int iW,int iT);

double DVImgGetTime(int iStream,int iZ,int iW,int iT);
double DVImgGetPhotoVal(int iStream,int iZ,int iW,int iT);
double DVImgGetMin(int iStream,int iZ,int iW,int iT);
double DVImgGetMax(int iStream,int iZ,int iW,int iT);
double DVImgGetMean(int iStream,int iZ,int iW,int iT);

double DVImgGetOriginX(int iStream);
double DVImgGetOriginY(int iStream);
double DVImgGetOriginZ(int iStream);

int    DVImgGetLensID(int iStream);

double DVImgGetWavelength(int iStream, int iWaveNum);

double DVImgGetIntenMin(int iStream, int iWaveNum);
double DVImgGetIntenMax(int iStream, int iWaveNum);
double DVImgGetIntenMean(int iStream, int iWaveNum);

double DVImgGetDisplayMin(int iStream, int iWaveNum);
double DVImgGetDisplayMax(int iStream, int iWaveNum);
double DVImgGetDisplayExp(int iStream, int iWaveNum);

int    DVImgGetNumCols(int iStream);    
int    DVImgGetNumRows(int iStream);    
int    DVImgGetNumZ(int iStream);
int    DVImgGetNumW(int iStream);
int    DVImgGetNumT(int iStream);

double DVImgGetLibVersion(void);
double DVImgGetBaseLibVersion(void);
int    DVImgGetBaseLibBuild(void);

char   *DVImgGetErrText(int iErrorNum);

double DVImgGetTimeStamp(int iStream);

// Header functions added 6/15/09 -KDC
int    DVImgGetFieldHdrSpace(int iStream);
int    DVImgGetNumFields(int iStream);
char   *DVImgGetFieldName(int iStream, int iFieldNum);
char   DVImgGetFieldType(int iStream, int iFieldNum);

int    DVImgSetFieldFloat(int iStream, char *cFieldName, float fValue);
int    DVImgSetFieldDouble(int iStream, char *cFieldName, double dValue);
int    DVImgSetFieldInt(int iStream, char *cFieldName, int iValue);
int    DVImgSetFieldChar(int iStream, char *cFieldName, char *cFieldContents);

float  DVImgGetFieldFloat(int iStream, char *cFieldName);
double DVImgGetFieldDouble(int iStream, char *cFieldName);
int    DVImgGetFieldInt(int iStream, char *cFieldName);
char   *DVImgGetFieldChar(int iStream, char *cFieldName);

// Image stream pointer positioning functions
int    DVImgMoveToZWT(int iStream, int iZ, int iW, int iT);

// Image IO functions
int    DVImgRead(int iStream, float *fArray);
int    DVImgWrite(int iStream, float *fArray);

// Auxiliary functions 
int    DVImgCheckFile(char *cFileName);
int    DVImgGetFileInfo(int iStream);

#ifdef WIN32
  // needed for external .def file used to export functions in the DLL
  #undef AFX_DATA
  #define AFX_DATA
#endif

#ifdef __cplusplus
}
#endif
