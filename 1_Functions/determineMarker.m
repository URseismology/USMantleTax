function marker = determineMarker(studyName, studyInterp)
    if strcmp(studyName, 'Liu')
        if strcmp(studyInterp, 'deepMLD')
            marker = 'h'; % Triangle for Liu deep
        else
            marker = 'p'; % Square for Liu shallow
        end
    elseif strcmp(studyName, 'Kreuger')
        marker = 'd'; % Pentagon for Kreuger
    elseif strcmp(studyName, 'Hua_PVG')
        marker = 's'; % Hexagon for Hua PVG
    elseif strcmp(studyName, 'Hua')
        marker = 's'; % Diamond for Hua
    else % Default case, for example, for Abt
        marker = '^'; % Circle for others
    end
end
