function ro = rollover(fig_hdl, varargin)
% Constructor for ROLLOVER objects

% Default = current figure
if nargin == 0
    fig_hdl = gcf;
end

% Copy constructor
if isa(fig_hdl, 'rollover')
    ro = fig_hdl;
    return
end

% Test validity of given figure handle
if ~ishandle(fig_hdl) || ~strcmp(get(fig_hdl, 'Type'), 'figure')
    error('Le premier argument doit etre un handle de figure valide !!');
end

% Default members
ro.Handles = [];
ro.StringsDefault = {};
ro.IconsDefault = {};
ro.StringsOver = {};
ro.IconsOver = {};
ro.CurrentButtonHdl = [];

% Read-only members
ro.Parent = fig_hdl; % Verrouillé à l'exécution

% Set figure's WindowButtonMotionFcn to activate rollover effect
set(fig_hdl, 'WindowButtonMotionFcn', 'roll');

% Declare class
ro = class(ro, 'rollover');

% Save application data into figure
setappdata(fig_hdl, 'rollover', ro);

% Parse property-value pairs
if nargin > 2
    ro = set(ro, varargin{:});
end