function year = get_publication_year(study)
    % This function returns the publication year for a given study
    switch study
        case 'Abt'
            year = '2010';
        case 'Hopper'
            year = '2018';
        case 'Hua'
            year = '2023';
        case 'Kreuger'
            year = '2021';
        otherwise
            year = '';
    end
end