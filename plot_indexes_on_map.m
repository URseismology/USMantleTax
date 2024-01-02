function plot_indexes_on_map(startIndex, endIndex, lats, lons, indices, uslatlim, uslonlim, colormapName)
    % Extract latitudes and longitudes for the specified range of data
    subsetLats = lats(indices(startIndex:endIndex));
    subsetLons = lons(indices(startIndex:endIndex));

    % Create a new figure
    figure(); clf;

    % Create a US map without state boundaries
    usamap(uslatlim, uslonlim);

    % Color code each station by its index in the specified range
    scatterm(subsetLats, subsetLons, 80, startIndex:endIndex, 'filled');

    % Use the specified colormap
    colormap(colormapName);
    colorbar;

    % Overlay the provincial borders
    geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'g', 'LineWidth', 1, 'FaceColor', 'none');

end
