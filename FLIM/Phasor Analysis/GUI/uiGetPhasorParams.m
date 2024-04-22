function params = uiGetPhasorParams
% Params struct keys
keys = {'ch', 'f', 'n', 'dt'};

% get values
prompts = {'Decay Channel(s)', 'Laser Frequency (MHz)', 'Harmonic', 'Temporal Bin Width (ns)', 'Ref. Lifetime (ns)'};
magOrder = {0, 6, 0, -9, -9};
defaults = {'1, 2', '80', '1', '0.0390625', '0'};
vals = inputdlg(prompts, 'Calibration Parameters', 1, defaults);

vals = cellfun(@(x,y) str2num(x)*10^y, vals', magOrder, 'UniformOutput', false);
args = [keys; vals(1:end-1)];
params = struct(args{:});
params.calibration.tau_ref = vals{end};

% Derive additional params
params.w = 2*params.n*pi*params.f;
params.calibration.m_ref = 1./sqrt(1+(params.w*params.calibration.tau_ref).^2);
params.calibration.phi_ref = atan(params.w*params.calibration.tau_ref);
