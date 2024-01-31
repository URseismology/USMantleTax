function [semblance] = compute_semblance(data, window_size, step)
% COMPUTE_SEMBLANCE Computes semblance for a seismic data matrix
%
% Inputs:
%   data: A matrix of seismic data
%   window_size: Size of the sliding window for semblance computation
%   step: Step size for sliding window
%
% Outputs:
%   semblance: Semblance of the data, with size [1, n_samples]

% Get the dimensions of the data
[n_traces, n_samples] = size(data);

% Initialize the semblance array
semblance = zeros(1, n_samples);

% Slide the window across the data
for t = 1:step:(n_samples - window_size + 1)
    % Extract the window from the data
    window_data = data(:, t:(t + window_size - 1));
    
    % Compute the semblance
    stack = mean(window_data, 1);
    stack_energy = sum(stack.^2);
    trace_energy = sum(window_data.^2, 2);
    semblance(t:t + window_size - 1) = semblance(t:t + window_size - 1) + (stack_energy / sum(trace_energy))^2;
end

% Normalize the semblance
semblance = semblance / max(semblance);

end