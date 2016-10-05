#!/bin/bash

#SBATCH --partition=super
#SBATCH --nodes=1
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mark.kittisopikul@utsouthwestern.edu

UUID=$1;
creationDir=$2;

module load matlab/$3

if [ -z "$USER"  ] || [ -z "$SLURM_JOB_ID" ]
then
    exit 111;
fi
if ! [ -d $creationDir ]
then
    exit 121;
fi

mkdir -p /tmp/$USER/$SLURM_JOB_ID

matlab -nodisplay -nosplash -noFigureWindows -logfile $creationDir/test.log  << MATLAB_CODE
try
    pc = parcluster('local');
    pc.JobStorageLocation = '/tmp/$USER/$SLURM_JOB_ID';
    //matlabpool(pc, getenv('SLURM_CPUS_ON_NODE'));
    load(['$creationDir' filesep 'in.mat']);
    out = cell(1,N);
    [out{:}] = fcn(input{:});
    save(['$creationDir' filesep 'out.mat'],'out');
catch(err)
    save(['$creationDir' filesep 'err.mat'],err);
end
exit;
MATLAB_CODE

if [ -n "$USER" ] && [ -n "$SLURM_JOB_ID" ]
then
    rm -rf /tmp/$USER/$SLURM_JOB_ID
fi
