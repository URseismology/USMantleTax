function plot_station_location(gridNumber, unsequencedData, sequencedData, uslatlim, uslonlim)
    % Find the index of the station in the unsequenced data
    unsequencedIndex = find(strcmp({unsequencedData.RFStruct.GridNumber}, gridNumber));
    
    % Find the index of the station in the reordered dataset
    reorderedIndex = find(strcmp(cellstr(sequencedData.grid_numbers_reordered), gridNumber));

    % Extract latitude and longitude of the station
    stationLat = unsequencedData.RFStruct(unsequencedIndex).Lat;
    stationLon = unsequencedData.RFStruct(unsequencedIndex).Lon;

    % Create a map and plot the location of the station
    figure(5); clf;
    usamap(uslatlim, uslonlim);
    geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'g', 'LineWidth', 1, 'FaceColor', 'none');
    geoshow(stationLat, stationLon, 'DisplayType', 'point', 'Marker', 'o', 'Color', 'r', 'MarkerSize', 10);
    
    title(['Station ', gridNumber, ': Unsequenced Index = ', num2str(unsequencedIndex), ', Reordered Index = ', num2str(reorderedIndex)]);
end
