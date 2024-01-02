function isInContourZone = check_contour_zone(matching_moho_depths, converted_times, contour_data_tPps, contour_data_tPss)
    % Check which points fall within the contour zones for tPps and tPss
    isInContourZone = false(size(matching_moho_depths));
    for i = 1:length(isInContourZone)
        for j = 1:length(contour_data_tPps)
            if inpolygon(matching_moho_depths(i), converted_times(i), contour_data_tPps{j}(1,:), contour_data_tPps{j}(2,:))
                isInContourZone(i) = true;
                break;
            end
        end
        if ~isInContourZone(i)  % Only check tPss if the point is not already in the tPps contour zone
            for j = 1:length(contour_data_tPss)
                if inpolygon(matching_moho_depths(i), converted_times(i), contour_data_tPss{j}(1,:), contour_data_tPss{j}(2,:))
                    isInContourZone(i) = true;
                    break;
                end
            end
        end
    end
end