% Add the current directory and subdirectories to the path
addpath(genpath(pwd));

%Load datasets
unsequencedData = load('RFDataWithCoords_raw_unfiltered.mat', 'RFStruct');
sequencedData = load('positives_full.mat');  % positives dataset

% Define time vector and data matrix for sequenced data
t = sequencedData.time_vector;
imageRFs_sequenced = sequencedData.dataset_imageRFs_pos_reordered;

%%
%Plot data with polarity correction
h1 = figure(1); clf;
ax1 = axes('Parent', h1); % Create axes in the figure
hold on;
for ii = 1:size(imageRFs_sequenced, 1)
    trace = imageRFs_sequenced(ii, :) - mean(imageRFs_sequenced(ii, :));
    trace_norm = trace / max(abs(trace));
    trace_norm = -trace_norm;

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;
    
    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), [1 1 1], 'none', 1,0.8);
    jbfill(t(positives), yvals(positives), zeroLine(positives), [0 0 1], 'none', 1,0.8);
    % plot(t, yvals, 'k');  % Commented out to remove the black zero line
end
xlim([6 30]);
ylim([-1, size(imageRFs_sequenced, 1) + 2]);  % Adjusted for spacing
xlabel('Time (s)');
ylabel('Ordered Station Index');
camroll(270);

% Ensure the plot is fully rendered
drawnow;
% Add a rectangle around the axes area
pos = get(ax1, 'Position');
annotation('rectangle', pos, 'EdgeColor', 'k', 'LineWidth', 1);

%%
%Extract and plot unsequenced data
numTraces = length(unsequencedData.RFStruct);
imageRFs_unsequenced = zeros(numTraces, length(t));
for ii = 1:numTraces
    imageRFs_unsequenced(ii, :) = unsequencedData.RFStruct(ii).Data;
end

h2 = figure(2); clf;
ax2 = axes('Parent', h2); % Create axes in the figure
hold on;
for ii = 1:numTraces
    trace_norm = imageRFs_unsequenced(ii, :) / max(abs(imageRFs_unsequenced(ii, :)));
    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    positives = trace_norm > 0;
    negatives = trace_norm < 0;
    
    jbfill(t(positives), yvals(positives), zeroLine(positives), [0 0 1], 'none', 1, 0.8);
    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), [1 0 0], 'none', 1, 0.8);
    % plot(t, yvals, 'k');  % Commented out to remove the black zero line
end
xlim([6 30]);
ylim([-1, numTraces + 2]);  % Adjusted for spacing
xlabel('Time (s)');
ylabel('Unsequenced Station Index');
camroll(270);

% Ensure the plot is fully rendered
drawnow;
% Add a rectangle around the axes area
pos = get(ax2, 'Position');
annotation('rectangle', pos, 'EdgeColor', 'k', 'LineWidth', 1);