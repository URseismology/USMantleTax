% Function to find and plot station by name
function plotStationLocationByName(stationName, stationNames, stationLats, stationLons, idx, controlPoints)
    % Find the index of the station
    stationIndex = find(strcmp(stationNames, strtrim(stationName)));
    
    if isempty(stationIndex)
        disp('Station name not found.');
        return;
    end
    
    % Get the cluster index for the station
    clusterIndex = idx(stationIndex);
    
    % Plot the map with the station location
    figure; usamap('conus'); hold on;
    
    % Add coastlines to the map
    states = shaperead('usastatelo', 'UseGeoCoords', true);
    geoshow([states.Lat], [states.Lon], 'DisplayType', 'polygon', 'FaceColor', [1 1 1])
    
    % Plot the station as a colored marker
    scatterm(stationLats(stationIndex), stationLons(stationIndex), 100, controlPoints(clusterIndex, :), 'filled');
    
    % Annotate the station name on the map
    textm(stationLats(stationIndex), stationLons(stationIndex), stationName, 'VerticalAlignment', 'bottom');
    
    % Set the title and other properties of the map
    title(sprintf('Location of Station %s (Cluster %d)', stationName, clusterIndex));
    
    hold off;
end
