% Load the dataset
data_table = readtable('merged_seismic_stations_data.csv');

% Convert 'Accepted' column to string if it's not already
if iscell(data_table.Accepted)
    data_table.Accepted = string(data_table.Accepted);
end

% Define the parameters for analysis
parameters = {'BazRatio', 'EpiDistRatio'};
accepted_logical = data_table.Accepted == "Yes";

% Initialize p-values array
p_values = zeros(1, length(parameters));

% Perform Mann-Whitney U test and plot histograms with KDE
figure(1);clf

for i = 1:length(parameters)
    subplot(length(parameters), 1, i);
    
    % Extracting data for the parameter
    accepted_data = data_table.(parameters{i})(accepted_logical);
    rejected_data = data_table.(parameters{i})(~accepted_logical);

    % Mann-Whitney U test
    [p, ~, ~] = ranksum(accepted_data, rejected_data);
    p_values(i) = p;

    % Histograms
    yyaxis left; % For histograms
    histogram(accepted_data, 'FaceColor', 'green', 'EdgeColor', 'black', 'Normalization', 'probability');
    hold on;
    histogram(rejected_data, 'FaceColor', 'red', 'EdgeColor', 'black', 'Normalization', 'probability');
    ylabel('Probability Density');
    set(gca, 'FontSize', 12); % Set font size for tick labels
    set(gca, 'YColor', 'black'); % Set y-axis label and tick colors to black

    % Drawing mean lines for histograms
    xline(mean(accepted_data), '--g', 'LineWidth', 2);
    xline(mean(rejected_data), '--r', 'LineWidth', 2);

    % KDE plots - excluding negative numbers and with continuous lines
    yyaxis right; % For KDE plots
    [f, xi] = ksdensity(accepted_data, 'Support', 'positive');
    plot(xi, f, 'LineWidth', 2, 'Color', 'green');
    [f, xi] = ksdensity(rejected_data, 'Support', 'positive');
    plot(xi, f, '-','LineWidth', 2, 'Color', 'red');
    set(gca, 'YColor', 'none'); % Make right y-axis ticks invisible

    % Set x-axis limits based on the parameter
    if strcmp(parameters{i}, 'EpiDistRatio')
        xlim([0 1]); % Set x-axis limit for EpiDistRatio
    elseif strcmp(parameters{i}, 'BazRatio')
        xlim([0 0.4]); % Set x-axis limit for BazRatio
    end

    % X-axis label
    xlabel(parameters{i});
    set(gca, 'FontSize', 17); % Set font size for tick labels

    % Title and legend
    title(['Histogram and KDE of ', parameters{i}, ' with Means']);
    legend(['Accepted'], ['Rejected'], 'Location', 'best');

    % Annotation for p-value
    annotation('textbox', [0.7, 0.8 - 0.4*(i-1), 0.1, 0.1], 'String', ['p-value: ' num2str(p_values(i), '%.1e')], 'FitBoxToText', 'on', 'BackgroundColor', 'white');

    hold off;
end

% Adjust layout
sgtitle('Histograms and KDEs of Key Parameters');
print(figure(1),'/Users/evets/Desktop/Earthscope/Figures/S7','-vector','-dpdf','-r0');



