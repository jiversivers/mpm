%% Imaging Parameters
% Load an IRF
[fn,fp] = uigetfile({'*.sdt'});
params = uiGetPhasorParams;
IRF = double(SDTtoTimeStack([fp filesep fn], 'channel', params.ch));
shift = find(sum(IRF, [1 2 3]) == max(sum(IRF, [1 2 3]), [], 4));

% Get IRF calib
params = IRFCalibration(IRF, params, 'off');

%% Simulate data
% Simulation parameters
I0 = 1; % peak photon counts
noisy = .025; % Max noise as percent of peak
blurry = 0.5;
decaySize = [1 1 1];

% Fluorophores
tau = [.2e-9 .4e-9];
alpha = [1 0];
alpha = alpha/(sum(alpha)); % normalize

% Simulate decay for each region
t = params.dt*(1/2:1:params.T-1/2);
t = shiftdim(t, 1-numel(decaySize)) .* ones(decaySize);
pure = zeros(decaySize);
for ii = 1:numel(alpha)
    pure = pure+(I0*alpha(ii)*exp(-t/tau(ii)));
end

% Add noise to decay (round to maintain discrete counts)
noised = simBlurryNoisyDecay(pure, 0, 0);

% Shift decay
really = cat(4, zeros([size(noised, [1 2 3]), shift]), noised(:,:, :, 1:end-shift));

% Convolve the two
really = convn(really, IRF);
really = really(:,:,:,shift:params.T+shift-1);

%% Test phasors
% Prepare figure and axes
% fig = figure;
% ax = axes;
% hold(ax, 'on');

%%
% Calculate phasor of pure decay
[g, s] = getPhasorCoords(pure, params, 'off');
% plotPhasor(ax, g, s);
% phasorView(ax, [], [], [], 'circleoff');
% addMultiLifetimeLine(tau, params.w);
% simOut.pure = getPhasorOutputs(g, s, params);
% %%
% % Calculate the phasor coords for the convovled
[g, s] = getPhasorCoords(really, params);
% plotPhasor(ax, g, s);
% phasorView(ax);
% addMultiLifetimeLine(tau, params.w);
% simOut.real = getPhasorOutputs(g, s, params);
% 
% % Compare to pure decay parameters
% simOut.groundTruth.Tau = sum(alpha.*tau);
[gb, sb] = lifetimetophasor(tau, params.w);
G = sum(alpha.*gb);
S = sum(alpha.*sb);

hold on
scatter(g, s, 'rx')
scatter(G, S)
viscircles(gca, [0.5, 0], 0.5);
yline(0); xline(0);

% Function that actually draws the values from the decay function
function pure = simDecayCurve(alpha, tau, t, I0, decaySize)
% alpha and tau are the species' parameters and must be the same length.
% t is a vector of timepoints to sample
% I0 is a scalar of max photon counts
% size is the size of the decay "image" to return
t = shiftdim(t, 1-numel(decaySize)) .* ones(decaySize);
pure = zeros(decaySize);
for ii = 1:numel(alpha)
    pure = pure+I0*alpha(ii)*exp(-t/tau(ii));
end
end

% Function that blurs and noises an input decay
function blurryNoisy = simBlurryNoisyDecay(pure, noisiness, blurriness)
% Pure is a decay curve that matches the rounded exponential decay purely
% Noisiness is a scalar, typically between 0 and 1 though not bound, that
% descrbies what inensity relative to the max photon count that noise is
% limited to.
% Blurriness is a scalar that describes the standard deviation of the
% gaussian filter used to blur the decay
blurryNoisy = pure+2*noisiness*max(pure, [], 'all')*(rand(size(pure))-0.5);
if blurriness > 0
    blurryNoisy = imgaussfilt(blurryNoisy, blurriness, 'Padding', 'symmetric');
end
blurryNoisy(blurryNoisy<0) = 0;
end



