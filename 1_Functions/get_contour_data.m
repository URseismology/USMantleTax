% Define the functions used in the script
function contour_data = get_contour_data(C)
    % Convert contour data to a format usable by inpolygon
    contour_data = {};
    start_idx = 1;
    while start_idx < size(C, 2)
        num_points = C(2, start_idx);
        contour_data{end+1} = C(:, start_idx+1:start_idx+num_points);  % Store x and y data for each contour line
        start_idx = start_idx + num_points + 1;
    end
end