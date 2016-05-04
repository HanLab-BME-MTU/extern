function [Xia,Yia,Uia,Via] = piv2InterpolateDisplacement(IAinfo,X,Y,U,V)
% piv2InterpolateDisplacement - returns displacements interpolated for IA grid
%
% Inputs: 
%    IAinfo ... structure containing information about IA (element of array created by IAinfo)
%    X,Y ... positions, for which displacement is known
%    U,V ... horizontal and vertical displacement
%
% Outputs:
%    Xia,Yia ... coordinaes of pixels inside IA
%    Uia,Via ... displacement interpoated for these pixels

%% 1. Get IA mesh
[Xia,Yia] = meshgrid(IAinfo.minX:IAinfo.maxX,IAinfo.minY:IAinfo.maxY);

%% 2. Interpolate displacement for pixels inside IA
Uia = interp2(X,Y,U,Xia,Yia,'*spline');
Via = interp2(X,Y,V,Xia,Yia,'*spline');