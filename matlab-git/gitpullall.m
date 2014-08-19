function gitpullall(basedir)
% gitpullall does a git pull in all subdirectories of basedir using .git
%
% If basedir is not specified, then the current directory is used.
% gitpullall is intended to be used in HOME/matlab
%
% Mark Kittisopikul
% August 19, 2014
% UT Southwestern
%
% license.txt of git.m does not apply to this file
%

    % save the currrent directory to return to when finished
    originaldir = pwd;

    if(nargin < 1)
        basedir = pwd;
    else
        % if a basedir is given, cd to it and obtain fullpath
        cd(basedir);
        basedir = pwd;
        disp(['* Changing to ' basedir ]);
    end

    D  = dir(basedir);
    % get directories only
    D = D([D.isdir]);
    
    for d = 1:length(D)
        % only go to directories with a .git folder
        if( exist([D(d).name filesep '.git'],'dir') )
            % cd to dir, pull, and then cd to base
            cd(D(d).name);
            disp(['    * Changing to ' D(d).name ' in ' basedir]);
            git('pull');
            disp(' ');
            cd(basedir);
        end
    end
    
    disp(['* Changing back to original directory: ' originaldir]);
    cd(originaldir);
    disp(' ');
    disp(['Done. If you do not see errors above, then you are up-to-date.']);

end
