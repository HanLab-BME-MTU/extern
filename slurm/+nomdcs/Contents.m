% The contents of this folder are used under exclusive license from
% THE MATHWORKS to UT Southwestern and are not to be redistributed.
%
% This contains code for interfacing Matlab with SLURM. For multinode jobs
% Matlab Distributed Computing Server is needed. For single node
% jobs a local license is sufficient. This code enables use of a single
% node without checking out MDCS licenses.
%
% Usage
% pc = parcluster('nucleus_r2015a');
% pc.CommunicatingSubmitFcn = @nomdcs.communicatingSubmitFcn;
% batch(pc,@() pause(3),0,'Pool',3)

