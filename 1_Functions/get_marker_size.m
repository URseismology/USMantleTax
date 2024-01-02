function size = get_marker_size(study)
    switch study
        case 'Hopper et al. 2018'
            size = 7;
        case 'Abt et al. 2010'
            size = 14;
        case 'Hua et al. 2023'
            size = 14;
        case 'Kreuger et al. 2021'
            size = 14;
        case 'Liu and Shearer 2021'
            size = 14;
        otherwise
            error('Unknown study size');
    end
end