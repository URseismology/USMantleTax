%  Reference: Carr, Olugboji(2024) 
%  A Taxonomy of Upper Mantle Stratification in the US.
%
%  Copyright (C) 2024, URSeismology lab
%  For more information: xxx
%  Author: Carr Steve
% 
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published
%  by the Free Software Foundation, either version 3 of the License, or
%  any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details: http://www.gnu.org/licenses/
%%
%First run a0_setupparameters
mntDir = '/gpfs/fs2/scratch/tolugboj_lab/';
outdir = [mntDir '/Prj7_RadonT/Prj7_US_Earthscope/2_Data/'];
addpath(outdir);

%list of seismic network and station names
network_station_pairs = {'CI','CWC'};

% loop over network and station pairs
for i = 1:size(network_station_pairs, 1)
    netname = network_station_pairs{i, 1};
    staname = network_station_pairs{i, 2};
    MTCfile = [RFDIR netname '_' staname '_RF.mat']; %or noCrust
    RFmat = load(MTCfile, 'time','radRF','qRF','svRF','bin');

    R = RFmat.radRF';
    t = RFmat.time;
    rayP = RFmat.bin(:, 2)';
    dist= RFmat.bin(:, 1);

    dt = t(2) - t(1);

    % cut data
    tBegin = find(t >  1, 1);
    tEnd   = find(t > 30, 1);
    t = t(tBegin:tEnd);
    R = R(:, tBegin:tEnd);
    
    %Remove bad data
    for k = 1:size(R, 1)

        % Find rows with any NaN values
        nanRows = any(isnan(R), 2);

        % Remove those rows
        R = R(~nanRows, :);

        % remove the corresponding elements from rayP as well
        rayP = rayP(~nanRows);
    end

    dq = (qmax-qmin)/(nq-1);  
    qs = qmin+dq*(0:1:nq-1);   
    dir_path = [localBaseDir 'Prj7_RadonT/Prj7_US_Earthscope/fista_results/' netname '_' staname '/'];

    % ensure the directory existsfactor = 2; 
    if ~exist(dir_path, 'dir')
        mkdir(dir_path);
    end

    if use_fista

        % parameters of the FISTA's algorithm
        lambda_vals = linspace(0.8, 2.7, 20);
        maxiter_vals = 20;

        tic 
       % try different values of the parameters
        results = cell(size(lambda_vals, 2), 3);
        parfor counter = 1:numel(lambda_vals)
            lambda = lambda_vals(counter);
            m_out = sparse_inverse_radon_fista(R', dt, rayP, qs, ...
                    fmin, fmax, lambda, mu, maxiter_vals);
            results(counter,:) = {lambda, maxiter_vals, m_out};
        end
        toc
        
        % save the results
        filename = [dir_path 'fista_on_' netname '_' staname '.mat'];
        saveResults(results, filename);
    end
end