function contaminated = isContaminated(depth, latitude, longitude, moho_lat, moho_long, tPps_avg, tPss_avg, Schultz_Pelkum_data)
    % Find the closest Moho point based on geographic location
    distances_moho = sqrt((latitude - moho_lat).^2 + (longitude - moho_long).^2);
    [~, closest_idx_moho] = min(distances_moho);
    
    % Find the closest Schultz and Pelkum data point based on geographic location
    distances_sp = sqrt((latitude - Schultz_Pelkum_data.latitude).^2 + (longitude - Schultz_Pelkum_data.longitude).^2);
    [~, closest_idx_sp] = min(distances_sp);
    
    % Retrieve Vp and Vs values from the Schultz and Pelkum data for the closest location
    Vp = Schultz_Pelkum_data.vp(closest_idx_sp);
    Vs = Schultz_Pelkum_data.vsv(closest_idx_sp);
    
    % Convert depth to time
    convert_depth_to_time = @(depth) depth .* ((1/Vs) - (1/Vp));
    converted_time = convert_depth_to_time(depth);
    
    % Check if the converted time falls within the range of tPps_avg or tPss_avg for the closest Moho point
    contaminated = (converted_time >= tPps_avg(closest_idx_moho) && converted_time <= tPss_avg(closest_idx_moho));
end
