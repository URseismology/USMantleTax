function marker = get_marker(study)
    switch study
        case 'Hopper et al. 2018'
            marker = '.';  % Circle
        case 'Abt et al. 2010'
            marker = '^';  % Triangle
        case 'Hua et al. 2023'
            marker = 's';  % Square
        case 'Kreuger et al. 2021'
            marker = 'd';  % Diamond
        case 'Liu and Shearer 2021'
            marker = 'p';
        otherwise
            error('Unknown study');
    end
end