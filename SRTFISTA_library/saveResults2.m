function saveResults2(dir_path, netname, staname, tBegin_val, tEnd_val, results)
    % Construct the filename
    filename = [dir_path 'fista_on_' netname '_' staname '_time_' num2str(tBegin_val) '_to_' num2str(tEnd_val) '.mat'];
    
    % Save results to the constructed filename
    save(filename, 'results');
end