function [  ] = ViewSection( Stream, Z, W, T )
%This is a very simple example script that takes a (Z, W, T) (0-indexed) coordinate and
%displays the section at that coordinate from within MATLAB.  It could
%potentially be used within another MATLAB script to troubleshoot a problem
%by displaying a section at a certain intermediary point in the operation.
%It assumes the library and image have already been loaded.

nx = DVImgGetNumCols(Stream);%get the dimensions
ny = DVImgGetNumRows(Stream);
Min = DVImgGetDisplayMin(Stream, W);
Max = DVImgGetDisplayMax(Stream, W);

Section = DVImgRead(Stream,Z,W,T);%read the specified section

figure('Units','pixels','Position',[100 100 nx ny])
imagesc(Section,[Min Max]);
set(gca,'Position',[0 0 1 1])