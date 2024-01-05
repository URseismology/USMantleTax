function matching_moho_depths = find_matching_moho_depths(data_lat, data_long, moho_lat, moho_long, moho_depth)
    matching_moho_depths = zeros(size(data_lat));
    for i = 1:length(data_lat)
        distances = sqrt((data_lat(i) - moho_lat).^2 + (data_long(i) - moho_long).^2);
        [~, closest_idx] = min(distances);
        matching_moho_depths(i) = moho_depth(closest_idx);
    end
end
