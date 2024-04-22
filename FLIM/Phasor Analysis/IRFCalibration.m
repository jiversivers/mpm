function params = IRFCalibration(IRF, params, saveopt, varargin)
%{
This fucntion will calculate the necessary calibration for phasor
parameters based on an input or selectede IRF file. The IRF file can be
input as a .tiff stack, a raw .sdt file, or a file path to either of these
file types. If the argument is left blank/empty, you will be prompted via a
file selection GUI to navigate to a file meeting either of  the
aformentioned file specifications.

PARAMS is a structure containing imaging parameters relevant for IRF and
calibration and phasor calcualtion generally. The params structrue will be
nested into the calibration structure with additional fields
updated/appended that are calculated from the base parameters that must be
input (either manually or through the GUI). These fields area s follows:
f                   The frequency of the laser in MHz
dt                  The temporal bin width in ns
n                   The desired harmonic; a positive integer value
calibration         A nested structure that contains the below fields for
                    all channels of the input IRF:
    tau_ref         The refernece lifetime in ns
    m_ref*          The reference modulation
    phi_ref*        The reference phase
* Optional fields. If they are not included, they will be added.

Appended and/or updated within this structure will be the following fields
as calculated from the inputs:
w                   The angular frequency harmonic of the laser (in rad/s)
T                   The number of collection bins

The PARAMS output is a structure containing calibration and imaging
parmeter information to be used in phasor coordinate calculations. It
contains Phi and M calibrations both as a map of pixel-wise calibrations
and as mean calibration across the entire input IRF. If the input IRF is a
single decay it is treated as a 1x1xT image, resulting in both calibration
styles being equivalent. The calibration information can be asccesed
through the style type (map or mean) style. For example,
params.calibration.Map.Phi will return the map of pixel-wise phi
calibrations. The params structure will have the updated params
structure (as described above) as well.

After calibration is complete, the udpated parameter structure may be saved
as a .mat file. This behvaior can be controlled programatically thorugh the
SAVEOPT positional argument (details below) or through a UI pop-up.
To save, input one:                               'on', 'save', true, 1
To be asked, leave the option blank or input:     'ask'
To skip saving, input one:                        'off', 'nosave', false, 0

Additionally, name-value pairs may be added to that wil be passed to the
getPhasorCoords function. Check there and SDTtoTimeStack for more
information on what options may be input.
%}

% Imaging parameters
if ~exist('params', 'var') || isempty(params)
    params = uiGetPhasorParams;
end 

% Get IRF File
if ~exist('IRF', 'var') || isempty(IRF)
    [IRF_file, IRF_path] = uigetfile({'*.sdt'; '*.tif'}, '\\deckard\BMEG\Rajaram-Lab\Ivers,Jesse\Codes\MPM_Processing\FLIM Code\', 'Select calibration file...');
    IRF = [IRF_path filesep IRF_file];
end

if ~isnumeric(IRF)
    [~, IRF_File, ~] = fileparts(IRF);
    % Load IRF Decay
    IRF = double(SDTtoTimeStack(IRF, 'channel', params.ch, varargin{:}));
else
    IRF_File = 'IRF_Calibration';
end

% Secondary (derivative) parameters for calculation of calibration
params.T = size(IRF, ndims(IRF));   % Bin number
params.w = 2*pi*params.f*params.n; % Angular frequency harmonic of laser
params.calibration.m_ref = 1/sqrt(1 + (params.w*params.calibration.tau_ref)^2); % Reference modulation
params.calibration.phi_ref = atan(params.w*params.calibration.tau_ref); % Reference phase

% Get calibration map for each channel
for ii = 1:numel(params.ch)
    % Get raw coordinates
    [g_irf, s_irf] = getPhasorCoords(IRF(:,:,params.ch(ii),:), params, 'off', varargin{:});
    
    % Convert to polar coordinates
    [phi_irf, m_irf] = cart2pol(g_irf, s_irf);

    % Calibration finals
    params.calibration.map(ii).M = params.calibration.m_ref./m_irf;
    params.calibration.map(ii).Phi = params.calibration.phi_ref-phi_irf;

    %% Bulk
    [g_irf, s_irf] = getPhasorCoords(mean(IRF(:,:,params.ch(ii),:), [1 2], 'omitnan'), params, 'off');

    % Convert to polar coordinates
    [phi_irf, m_irf] = cart2pol(g_irf, s_irf);

    % Calibration finals
    params.calibration.mean(ii).M = params.calibration.m_ref/m_irf;
    params.calibration.mean(ii).Phi = params.calibration.phi_ref-phi_irf;
end

% Save Calibration
if~exist("saveopt", "var")
    saveopt = 'ask';
end
switch lower(saveopt)
    case {'save', 'on', true, 1}
        [sn, sp] = uiputfile('*.mat', 'Save new calibration...', [IRF_File, '.mat']); 
        save([sp sn], 'params', '-mat')
    case {'nosave', 'off', false, 0}
        return
    otherwise
        saveOpt = questdlg('Save new IRF Calibration?', 'Save...', 'Yes', 'No', 'Yes');
        if strcmp(saveOpt, 'Yes')
            [sn, sp] = uiputfile('*.mat', 'Save new calibration...', [IRF_File, '.mat']); 
            save([sp sn], 'params', '-mat')
        end
end