% Function to get color based on UMD
function color = get_color(UMD)
    switch UMD
        case 'UMD1'
            color = 'b';  % Blue for UMD1
        case 'UMD2'
            color = 'm';  % Red for UMD2
        case 'UMD3'
            color = 'g';  % Green for UMD3
        otherwise
            error('Unknown UMD');
            %color = 'b';  % Green for UMD3
    end
end