function phasorOutputs = getPhasorOutputs(G, S, w)

% Omega input can be calibration structure, parameter structure, or simply
% the value (in which case it does not meet this condition and is passed
% directly)
if isstruct(w)
    if any(contains(fields(w), 'w'))
        % For params struct
        w = w.w;
    else
        try
            % Calibration struct
            w=w.params.w;
        catch
            warning('No angular frequency found in parameter structure\nCalculating from frequency...\n')
            w=2*pi*w.params.f;
        end
    end
end

phasorOutputs.Phi = atan(S./G);
phasorOutputs.M = sqrt(G.^2 + S.^2);
phasorOutputs.Tau.Phi = (1/w)*tan(phasorOutputs.Phi);
phasorOutputs.Tau.M = (1/w)*sqrt(phasorOutputs.M.^(-2) - 1);