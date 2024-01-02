function nextCluster = determineNextCluster(currentCluster)
    maxCluster = 4; % Assuming the maximum cluster number is 4
    if currentCluster < maxCluster
        nextCluster = currentCluster + 1;
    else
        nextCluster = currentCluster; % If already at max, stay at max
    end
end
