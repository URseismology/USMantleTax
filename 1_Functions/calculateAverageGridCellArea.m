function avg_area = calculateAverageGridCellArea(lat, long, num_lat_bins, num_long_bins)
    % Define latitude and longitude bins
    lat_bins = linspace(min(lat), max(lat), num_lat_bins+1);
    long_bins = linspace(min(long), max(long), num_long_bins+1);
    
    % Initialize total area
    total_area = 0;
    
    % Loop through each grid cell and calculate its area
    for i = 1:num_lat_bins
        for j = 1:num_long_bins
            % Extract the bounds of the current grid cell
            lat1 = lat_bins(i);
            lat2 = lat_bins(i+1);
            lon1 = long_bins(j);
            lon2 = long_bins(j+1);
            
            % Calculate and accumulate the area of the current grid cell
            total_area = total_area + calculateGridCellArea(lat1, lat2, lon1, lon2);
        end
    end
    
    % Calculate the average area of the grid cells
    avg_area = total_area / (num_lat_bins * num_long_bins);
end

% Provided function to calculate the area of a grid cell
function area = calculateGridCellArea(lat1, lat2, lon1, lon2)
    R = 6371;  % Earth's radius in km
    lat1 = deg2rad(lat1);
    lat2 = deg2rad(lat2);
    lon1 = deg2rad(lon1);
    lon2 = deg2rad(lon2);
    
    vertical_distance = 2 * R * asin(sqrt(sin((lat2-lat1)/2)^2));
    horizontal_distance = 2 * R * asin(sqrt(cos(lat1)*cos(lat2)*sin((lon2-lon1)/2)^2));
    
    area = vertical_distance * horizontal_distance;
end