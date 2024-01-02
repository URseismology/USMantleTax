% Add the current directory and subdirectories to the path
addpath(genpath(pwd));

%Load datasets
load('full.mat');  % or negatives_full.mat
unsequencedData = load('RFDataWithCoords_raw_unfiltered.mat', 'RFStruct');
sequencedData = load('full.mat');  % positives dataset

% load('full.mat');  % or negatives_full.mat
% unsequencedData = load('RFDataWithCoords_negs_new.mat', 'RFStruct');
% sequencedData = load('negatives_full.mat');  %negatives dataset

% Define time vector and data matrix for sequenced and unsequenced data
t = time_vector;
imageRFs_sequenced = dataset_imageRFs_pos_reordered;

%%
%Plot data with polarity correction
h1 = figure(1); clf;
hold on;
for ii = 1:size(imageRFs_sequenced, 1)
    trace = imageRFs_sequenced(ii, :) - mean(imageRFs_sequenced(ii, :));
    trace_norm = trace / max(abs(trace));
      
    trace_norm = -trace_norm;

    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    negatives = trace_norm < 0;
    positives = trace_norm > 0;
    
    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), [1 0 0], 'none', 0);
    jbfill(t(positives), yvals(positives), zeroLine(positives), [0 0 1], 'none', 0);
    % plot(t, yvals, 'k');  % Commented out to remove the black zero line
end
xlim([6 30]);
ylim([-1, size(imageRFs_sequenced, 1) + 1]);  % Adjusted for spacing
xlabel('Time (s)');
ylabel('Ordered Station Index');
camroll(270);
%print(figure(1),'./Figures/figure4d','-vector','-dpdf','-r0');
%%
%Extract and plot unsequenced data
numTraces = length(unsequencedData.RFStruct);
imageRFs_unsequenced = zeros(numTraces, length(t));
for ii = 1:numTraces
    imageRFs_unsequenced(ii, :) = unsequencedData.RFStruct(ii).Data;
end

h2 = figure(2); clf;
hold on;
for ii = 1:numTraces
    trace_norm = imageRFs_unsequenced(ii, :) / max(abs(imageRFs_unsequenced(ii, :)));
    yvals = trace_norm + ii;
    zeroLine = ii * ones(size(t));
    positives = trace_norm > 0;
    negatives = trace_norm < 0;
    
    jbfill(t(positives), yvals(positives), zeroLine(positives), [0 0 1], 'none', 0);
    jbfill(t(negatives), yvals(negatives), zeroLine(negatives), [1 0 0], 'none', 0);
    % plot(t, yvals, 'k');  % Commented out to remove the black zero line
end
xlim([6 30]);
ylim([-1, numTraces + 1]);  % Adjusted for spacing
xlabel('Time (s)');
ylabel('Unsequenced Station Index');
camroll(270);
%print(figure(2),'./Figures/figure4d','-vector','-dpdf','-r0');
