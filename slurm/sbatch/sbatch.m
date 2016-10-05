function [ outStruct ] = sbatch( fcn, N, input, varargin )
%SBATCH Submit a batch job using slurm
%
% INPUT
% fcn   - function to execute remotely
% N     - number of outputs
% input - cell array of arguments to fcn
%
% OUTPUT
% outStruct  - containing the following fields
% .getStatus - function to retrieve the status
% .fetchOutput - function to load the output of the job
% .delete - function to delete job directory saved to disk
% .slurmStatus - function to query the status using sacct
% .diary - function to retrieve the standard output of the job
% .errors - function to retrieve the standard error of the job
% .getError - function to retrieve the MATLAB error that occurred
% .submissionOutput - string containing the output of sbatch on submission
% .uuid - universally unique identifier
% .creationDir - directory where job is stored
%
% FILES CREATED
% creationDir/in.mat - stores input to fcn
% creationDir/out.mat - stores output from fcn
% creationDir/job.mat - stores outStruct
% creationDir/out.log - stores standard output
% creationDir/test.log - unknown, created by SLURM?
%
% AUTHOR
% Mark Kittisopikul September 29th, 2016
% Jaqaman Lab
% UT Southwestern

if(nargin < 3)
    input = {};
end

% Create universally unique id via Java
uuid = char(java.util.UUID.randomUUID().toString());
uuid = ['sbatchMatlab_' uuid];
% Make sure uuid is not empty or this could go very badly
assert(~isempty(uuid),'sbatch:emptyUUID','UUID is empty');
creationDir = [pwd filesep uuid];
mkdir(uuid);
save([creationDir filesep 'in.mat'],'fcn','N','input');

% Determine shell script location should be sbatchMatlab.sh in the same
% folder as this file
shellScript = mfilename('fullpath');
shellScript = [shellScript 'Matlab.sh'];
% Submit job to SLURM
cmdin = [ ...
    'sbatch ' ...
    '--job-name="' func2str(fcn) '" ' ...
    '--output="' creationDir filesep 'out.log" ' ...
    ' --error="' creationDir filesep 'err.log" ' ...
    shellScript ' ' ...
    uuid ' '...  % First argument = UUID
    '"' creationDir '" ' ... % Second argument = creationDir
    version('-release') ... % Third argument = release (e.g. 2015a)
    ];
[outStruct.submissionStatus, cmdout] = system(cmdin);

% Create outStruct
outStruct.getStatus = @getStatus;
outStruct.fetchOutput = @fetchOutput;
outStruct.delete = @delete;
outStruct.slurmStatus = @slurmStatus;
outStruct.diary = @diary;
outStruct.errors = @errors;
outStruct.getError = @getError;
outStruct.submissionOutput = strtrim(cmdout);
outStruct.uuid = uuid;
outStruct.slurmID = regexp(outStruct.submissionOutput,'[0-9]*$','match');
outStruct.slurmID = outStruct.slurmID{end};
outStruct.creationDir = creationDir;
save([creationDir filesep 'job.mat'],'outStruct');

%% Accessory functions accessible via outStruct

    function status = getStatus()
        % GETSTATUS
        outExists = exist([creationDir filesep 'out.mat'],'file');
        if(outExists)
            status = 'finished';
        else
            [returnStatus,status] = system(['sacct -nPj' outStruct.slurmID ' --format=state']);
            if(returnStatus)
                % command failed
                status = 'not finished';
            else
                % command good
                status = strsplit(status);
                status = status{1};
                status = lower(status);
            end
        end
        if(nargout == 0)
            slurmStatus();
        end
    end
    function out = fetchOutput()
        % FETCHOUTPUT fetch sbatch output
        out = load([creationDir filesep 'out.mat']);
        out = out.out;
    end
    function delete()
        % DELETE sbatch files
        
        % Delete mat files
%         matFiles = dir([creationDir filesep '*.mat']);
%         matFiles = {matFiles.name};
        matFiles = strcat({'in','out','job'},'.mat');
        if(~isempty(matFiles))
            matFiles = strcat([creationDir filesep],matFiles);
            builtin('delete',matFiles{:});
        end
        
        % Delete log files
%        logFiles = dir([creationDir filesep '*.log']);
%        logFiles = {logFiles.name};
        logFiles = strcat({'test','out','err'},'.log');
        if(~isempty(logFiles))
            logFiles = strcat([creationDir filesep],logFiles);
            builtin('delete',logFiles{:});
        end
        
        % Remove Directory
        if(exist(creationDir,'dir'))
            rmdir([creationDir]);
        end
    end
    function slurmStatus()
        system(['sacct -j' outStruct.slurmID],'-echo');
    end
    function diary()
        type([creationDir filesep 'out.log']);
    end
    function errors()
        type([creationDir filesep 'err.log']);
    end
    function err = getError()
        errFile = [creationDir filesep 'err.mat'];
        if(exist(errFile,'file'))
            err = load([creationDir filesep 'err.mat']);
            err = err.err;
        else
            err = [];
        end
    end

end