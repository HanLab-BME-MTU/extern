function roll
% Implementation of the rollover feature via the figure's event
% WindowButtonMotionFcn

% Handle of the figure associated with WindowButtonMotionFcn event
if strcmp(get(gcbo, 'Type'), 'figure')
    fig_hdl = gcbo;
else
    error('The rollover feature must be implemented via a figure''s WindowButtonMotionFcn event !!');
end

% Mouse pointer over which control ?
current_object = hittest;

% Retrieve ROLLOVER object from current figure
ro = getappdata(fig_hdl, 'rollover');

% Test whether current_object is a button or not
if strcmp(get(current_object, 'Type'), 'uicontrol') && strcmp(get(current_object, 'Style'), 'pushbutton')
    % List of rollover-capable pushbuttons
    allowed_buttons = get(ro, 'Handles');

    % current_object belongs to previous list ?
    if isempty(find(current_object == allowed_buttons, 1))
        return
    end

    % Change ro.CurrentButtonHdl to current_object (and update label and
    % icon)
    set(ro, 'CurrentButtonHdl', current_object);
else
    % current_object is not a button
    % If ro.CurrentButtonHdl is not empty, pushbutton is being left by the
    % mouse pointer -> revert to default label and icon with
    % set(ro, 'CurrentButtonHdl', []);
    % Otherwise, do nothing (labels and icons of other buttons have already
    % been reverted back to normal)
    if ~isempty(get(ro, 'CurrentButtonHdl'))
        set(ro, 'CurrentButtonHdl', []);
    end
end