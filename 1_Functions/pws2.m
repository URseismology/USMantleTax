function [pwstack, mean_stack,phase_wts,wtdata] = pws2(data, p)
% PHASE_WEIGHTED_STACK Computes the phase-weighted stack of a seismic data matrix
%
% Inputs:
%   data: A matrix of seismic data 
%   p: The exponent used to weight the phase coherence.
%
% Outputs:
%   pwstack: The phase-weighted stack of the data, with size [1, n_samples].
%   phase_wts: The phase weights (coherence) applied to the data, with size [1, n_samples].
% Equations from Schimmel, M., & Paulssen, H. (1997). Noise reduction and detection of weak,coherent signals through phase-weighted stacks. 

n_params = size(data, 1);

% Compute the average trace as the reference signal
mean_stack = mean(data, 1);
mean_stack = mean_stack ./ max(mean_stack);

c_phase = zeros(n_params, length(mean_stack));

% Loop over all the ray parameters
for ii = 1:n_params
    % Compute the complex trace of the current parameter
    complex_trace = hilbert(data(ii, :));
    
    % Compute the instantaneous phase of the complex trace
    phase = angle(complex_trace);

    c_phase(ii, :) = exp(1i * phase); % complex Instantaneous phase for every rf trace
end

% Compute phase weights (coherency).
phase_wts = abs(mean(c_phase, 1)).^p; % Equation 2 in Schimmel and Paulssen (1997) coherency measure is independent of the amplitude

% Apply phase weighting to data
wtdata = data .* phase_wts; % Apply the absolute value only to the phase weights

% Stack data
pwstack = mean(wtdata);
% pwstack = pwstack / max(abs(pwstack));

end