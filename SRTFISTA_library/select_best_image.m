function [best_m_out,best_lambda,best_score] = select_best_image(results)
% Define the regions of interest in time (seconds)
regions_of_interest = {[6, 15], [12, 15], [23, 33], [30, 37], [38, 45]};

% Initialize arrays to store the energy metrics, sparsity values and lambda values for each image
energy_metrics = zeros(1, 10);
sparsity_values = zeros(1, 10);
lambda_values = zeros(1, 10);

% Calculate energy metrics, sparsity values and lambda values for each image
for i=1:6
    m_out = cell2mat(results(i,3));
    lambda = cell2mat(results(i,1));

    % Compute energy concentration in regions of interest
    %     energy_metric = 0;
    %     for region = regions_of_interest
    %         start_index = round((region{1}(1) - 2) * (size(m_out, 1) / (50 - 2)));
    %         end_index = round((region{1}(2) - 2) * (size(m_out, 1) / (50 - 2)));
    %         region_indices = start_index:end_index;
    %         region_energy = sum(sum(m_out(region_indices, :).^2));
    %         total_energy = sum(sum(m_out.^2));
    %         energy_metric = energy_metric + (region_energy / total_energy);
    %     end

    energy_metric = sum(sum(m_out.^2)) / sum(sum(m_out.^2));

    % Compute sparsity (L1 norm)
    sparsity = sum(abs(m_out), 'all');

    % Store the energy metric, sparsity value and lambda value
    energy_metrics(i) = energy_metric;
    sparsity_values(i) = sparsity;
    lambda_values(i) = lambda;
end

% Normalize the energy metrics, sparsity values and lambda values to the range [0,1]
energy_metrics_norm = (energy_metrics - min(energy_metrics)) / (max(energy_metrics) - min(energy_metrics));
sparsity_values_norm = (sparsity_values - min(sparsity_values)) / (max(sparsity_values) - min(sparsity_values));
lambda_values_norm = (lambda_values - min(lambda_values)) / (max(lambda_values) - min(lambda_values));

% Define the weights for the criteria (adjust these according to your needs)
w_energy = 0.6;  % weight for the energy concentration in the regions of interest
w_sparsity = 0.3;  % weight for the sparsity of the image
w_lambda = 0.1;  % weight for the lambda value (we want this to be low)

% Initialize the best score and corresponding best index
best_score = -Inf;
best_index = 0;

% Calculate scores for each image
for i=1:6
    % Compute the score
    score = w_energy * energy_metrics_norm(i) - w_sparsity * sparsity_values_norm(i) - w_lambda * lambda_values_norm(i);

    % If score is better than current best, update the best score and best index
    if score > best_score
        best_score = score;
        best_index = i;
    end
end

% Display the best index
disp(['Best image index: ', num2str(best_index)]);
% After getting best_index
best_m_out = cell2mat(results(best_index,3));
best_lambda = cell2mat(results(best_index,1));

% Calculate scores
scores = w_energy * energy_metrics_norm - w_sparsity * sparsity_values_norm - w_lambda * lambda_values_norm;

best_score = scores(best_index);

% Return the best index and best image
return;
end