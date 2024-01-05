function getStationDataByCluster(idx, stationNames, stationLats, stationLons, clusterNumber, filename)
    % This function gets the names, latitudes, and longitudes of stations that belong to the specified cluster
    % and saves them to a text file.
    %
    % Parameters:
    % idx - A vector of cluster indices for each station
    % stationNames - A cell array of station names
    % stationLats - A vector of station latitudes
    % stationLons - A vector of station longitudes
    % clusterNumber - The cluster number you're interested in
    % filename - The name of the text file to save the data
    
    % Find the indices of the stations that belong to the cluster
    clusterIndices = find(idx == clusterNumber);
    
    % Get the data for the stations in that cluster
    stationNamesInCluster = stationNames(clusterIndices);
    stationLatsInCluster = stationLats(clusterIndices);
    stationLonsInCluster = stationLons(clusterIndices);
    
    % Open the file for writing
    fileID = fopen(filename, 'w');
    
    % Check if the file was opened successfully
    if fileID == -1
        error('Unable to open file %s for writing.', filename);
    end
    
    % Write the data to the file
    for i = 1:length(clusterIndices)
        fprintf(fileID, '%s\t%f\t%f\n', stationNamesInCluster{i}, stationLatsInCluster(i), stationLonsInCluster(i));
    end
    
    % Close the file
    fclose(fileID);
end

