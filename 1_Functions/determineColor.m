function color = determineColor(studyName, studyInterp)
    if strcmp(studyName, 'Hua_PVG') || (strcmp(studyName, 'Kreuger') && strcmp(studyInterp, 'PVG'))
        color = 'b'; % Blue for PVG
    else
        color = 'r'; % Red for NVG and others
    end
end
