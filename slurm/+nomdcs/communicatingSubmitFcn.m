function communicatingSubmitFcn(cluster, job, props)
%COMMUNICATINGSUBMITFCN Submit a communicating MATLAB job to a SLURM cluster
%
% Set your cluster's CommunicatingSubmitFcn to this function using the following
% command:
%     set(cluster, 'CommunicatingSubmitFcn', @communicatingSubmitFcn);
%
% See also parallel.cluster.generic.communicatingDecodeFcn.
%

% Copyright 2010-2012 The MathWorks, Inc.

% Store the current filename for the errors, warnings and dctSchedulerMessages
currFilename = mfilename;
if ~isa(cluster, 'parallel.Cluster')
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The function %s is for use with clusters created using the parcluster command.', currFilename)
end

decodeFunction = 'parallel.cluster.generic.communicatingDecodeFcn';

if ~cluster.HasSharedFilesystem
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The submit function %s is for use with shared filesystems.', currFilename)
end

if ~strcmpi(cluster.OperatingSystem, 'unix')
    error('parallelexamples:GenericSLURM:SubmitFcnError', ...
        'The submit function %s only supports clusters with unix OS.', currFilename)
end

% The job specific environment variables
% Remove leading and trailing whitespace from the MATLAB arguments
matlabArguments = strtrim(props.MatlabArguments);
variables = {'MDCE_DECODE_FUNCTION', decodeFunction; ...
    'MDCE_STORAGE_CONSTRUCTOR', props.StorageConstructor; ...
    'MDCE_JOB_LOCATION', props.JobLocation; ...
    'MDCE_MATLAB_EXE', props.MatlabExecutable; ...
    'MDCE_MATLAB_ARGS', matlabArguments; ...
    'MDCE_DEBUG', 'true'; ...
    'MLM_WEB_LICENSE', props.UseMathworksHostedLicensing; ...
    'MLM_WEB_USER_CRED', props.UserToken; ...
    'MLM_WEB_ID', props.LicenseWebID; ...
    'MDCE_LICENSE_NUMBER', props.LicenseNumber; ...
    'MDCE_STORAGE_LOCATION', props.StorageLocation; ...
    'MDCE_CMR', cluster.ClusterMatlabRoot; ...
    'MDCE_TOTAL_TASKS', num2str(props.NumberOfTasks)};
% Trim the environment variables of empty values.
nonEmptyValues = cellfun(@(x) ~isempty(strtrim(x)), variables(:,2));
variables = variables(nonEmptyValues, :);
% Set the remaining non-empty environment variables
for ii = 1:size(variables, 1)
    setenv(variables{ii,1}, variables{ii,2});
end

% The script name is communicatingJobWrapper.sh
scriptName = 'communicatingJobWrapper.sh';
% The wrapper script is in the same directory as this file
dirpart = fileparts(mfilename('fullpath'));
scriptName = fullfile(dirpart, scriptName);

% Choose a file for the output. Please note that currently, JobStorageLocation refers
% to a directory on disk, but this may change in the future.
logFile = cluster.getLogLocation(job);

jobName = sprintf('Job%d', job.ID);
% SLURM jobs names must not exceed 15 characters
maxJobNameLength = 15;
if length(jobName) > maxJobNameLength
    jobName = jobName(1:maxJobNameLength);
end
dctSchedulerMessage(5, '%s: Generating command for task %i', currFilename, ii);
commandToRun = getSubmitString(jobName, logFile, scriptName, ...
    props);

% Now ask the cluster to run the submission command
dctSchedulerMessage(4, '%s: Submitting job using command:\n\t%s', currFilename, commandToRun);
try
    % Make the shelled out call to run the command.
    [cmdFailed, cmdOut] = system(commandToRun);
catch err
    cmdFailed = true;
    cmdOut = err.message;
end
if cmdFailed
    error('parallelexamples:GenericSLURM:SubmissionFailed', ...
        'Submit failed with the following message:\n%s', cmdOut);
end

dctSchedulerMessage(1, '%s: Job output will be written to: %s\nSubmission output: %s\n', currFilename, logFile, cmdOut);

jobIDs = extractJobId(cmdOut);
% jobIDs must be a cell array
if isempty(jobIDs)
    warning('parallelexamples:GenericSLURM:FailedToParseSubmissionOutput', ...
        'Failed to parse the job identifier from the submission output: "%s"', ...
        cmdOut);
end
if ~iscell(jobIDs)
    jobIDs = {jobIDs};
end

% set the job ID on the job cluster data
cluster.setJobClusterData(job, struct('ClusterJobIDs', {jobIDs}));

function submitString = getSubmitString(jobName, logFile, command, ...
    props)
%GETSUBMITSTRING Gets the correct sbatch command for a SLURM cluster

% Copyright 2010-2012 The MathWorks, Inc.

additionalSubmitArgs = getAdditionalSubmitArguments(props)

submitString = sprintf('sbatch --job-name=%s --output="%s" %s "%s"', ...
       jobName, logFile, additionalSubmitArgs, command);


function asa = getAdditionalSubmitArguments(props)

% Copyright 2010-2012 The MathWorks, Inc.

asa = '';


%% REQUIRED

asa = [asa ' --ntasks=' num2str(props.NumberOfTasks)];


%% OPTIONAL

% EMail notification
ea = ClusterInfo.getEmailAddress();
if isempty(ea)==false
    asa = [asa ' --mail-type=ALL --mail-user=' ea];
end

% GPU
gpu = ClusterInfo.getUseGpu();
if gpu==true
    asa = [asa ' -p GPU'];
else
   % Check for partition
    qn = ClusterInfo.getQueueName();
    if isempty(qn)==false
        asa = [asa ' -p ' qn];
   end
end

% Walltime
wt = ClusterInfo.getWallTime();
if isempty(wt)==false
    asa = [asa ' -t ' wt];
end

% NNode
nnodes = ClusterInfo.getNNode();
useMDCS = false;
if ~isempty(nnodes) && nnodes ~= 1
    asa = [asa ' -N ' num2str(nnodes)];
    useMDCS = true;
else
    % mkitti: 2016_10_04, restrict to one node only
    asa = [asa ' -N 1'];
end

% mkitti: 2016_10_04, no mdcs license needed if nnodes == 1
if(useMDCS)
    [~,flag] = strtok(props.MatlabArguments,'p');
    if strcmp(flag,'parallel')
        w = props.NumberOfTasks;
    else
        w = 1;
    end


    % Specification of MDCS licenses which must be allocated to this
    % job.  The /etc/slurm/slurm.conf file must list
    %
    %   # MDCS licenses
    %   Licenses=mdcs:128
    %
    asa = [asa ' --licenses=mdcs*' num2str(w)];
end

% Catch-all
udo = ClusterInfo.getUserDefinedOptions();
if isempty(udo)==false
    asa = [asa ' ' udo];
end

if isempty(asa)==false
    asa = strtrim(asa);
end

