function F = set(varargin)


% I AM SO SORRY FOR MESSING WITH
% INTERNAL FUNCTIONS; BUT I REALLY
% WANT TO BE ABLE TO DO set([]) and set;

% Suppressing warning https://git.biohpc.swmed.edu/danuser/extern/commit/fdcc3e912cd4dca3f8fd70f08329db12315a40e7#note_4957
%warning('yalmip:double:set','Yalmip in extra intercepted set!');

% CATCH MY CASES
switch nargin
    case 0
        F = lmi;
        return
    case 1
        if isempty(varargin{1})
            F = lmi;
            return
        end
    otherwise
end

% ... Leave for standard code
try
    if nargout == 0
        builtin('set',varargin{:});
    else
        F = builtin('set',varargin{:});
    end
catch
    error(lasterr);
end

