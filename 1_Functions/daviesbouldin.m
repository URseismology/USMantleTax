function dbIndex = daviesbouldin(X, labels)
% X is the data matrix (n x m, with n samples and m features)
% labels are the cluster labels for each sample

k = max(labels); % number of clusters
cluster_k = cell(k, 1);
centroids = zeros(k, size(X, 2));

% Compute cluster centroids and dispersion (S)
S = zeros(k, 1);
for i = 1:k
    cluster_k{i} = X(labels == i, :);
    centroids(i, :) = mean(cluster_k{i}, 1);
    distances = pdist2(cluster_k{i}, centroids(i, :));
    S(i) = mean(distances);
end

% Compute the similarity measure (R) between each pair of clusters
R = zeros(k, k);
for i = 1:k
    for j = i+1:k
        Mij = pdist2(centroids(i, :), centroids(j, :));
        R(i, j) = (S(i) + S(j)) / Mij;
        R(j, i) = R(i, j); % The matrix is symmetric
    end
    R(i, i) = 0; % The diagonal should be 0 as similarity with itself is not defined
end

% Compute Davies-Bouldin index
Ri = max(R, [], 2);
dbIndex = mean(Ri(Ri ~= Inf)); % Avoid Inf values which occur if Mij is 0

% Return NaN if there's only one cluster (k=1) since DB index isn't defined for a single cluster
if k == 1
    dbIndex = NaN;
else
    dbIndex = mean(Ri);
end
end




