% This folder contains a script to run a function on a remote compute node using Slurm. There are two parts:
% 1) A MATLAB function sbatch.m that saves the input arguments and parameters similar to the batch command.
% 2) A BASH script for Slurm job submission. This is the script that actually calls MATLAB and passes in the appropriate commands
%    to load the input, run the function, and save the output
