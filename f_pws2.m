function [pwstack_semblance] = f_pws2(data, p, window_size, step)
% f_pws2 Computes the phase-weighted stack incorporating semblance
%
% Inputs:
%   data: A matrix of seismic data
%   p: The exponent used to weight the phase coherence.
%   window_size: Size of the sliding window for semblance computation
%   step: Step size for sliding window
%
% Outputs:
%   pwstack_semblance: The phase-weighted stack of the data incorporating semblance, with size [1, n_samples].

    % Compute the phase-weighted stack
    [~, ~, phase_wts, wtdata] = pws2(data, p);

    % Compute the semblance
    semblance = compute_semblance(data, window_size, step);

    % Combine semblance with phase weights
    combined_weights = phase_wts .* semblance;

    % Apply combined weights to data
    wtdata_semblance = data .* combined_weights;

    % Stack data with semblance weights
    pwstack_semblance = mean(wtdata_semblance);
%     pwstack_semblance = pwstack_semblance / max(abs(pwstack_semblance));
    pwstack_semblance = pwstack_semblance';

end