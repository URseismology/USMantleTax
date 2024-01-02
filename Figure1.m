% Load NVG depths from Excel file
addpath(genpath(pwd));

filename1 = 'Hopper2018.xlsx';
filename2 = 'Abt_NVG_2010.txt';
filename3 = 'HUA_MLD.txt';
filename4 = 'HUA_PVG.txt';
filename5 = 'Liu_Shallow_MLD.txt';
filename6 = 'Liu_LAB.txt';
filename7 = 'Liu_deep_MLD.txt';

% Read UMD info from Hopper data
HopperMLD = xlsread(filename1); 
latitude_Hopper_MLD = HopperMLD(:,1);
longitude_Hopper_MLD = HopperMLD(:,2);
depth_Hopper_MLD = HopperMLD(:,3);

% Read UMD info from Abt data
AbtMLD = readtable(filename2, 'Delimiter', '\t');
latitude_Abt = AbtMLD.Lat;
longitude_Abt = AbtMLD.Long;
depth_Abt = AbtMLD.Depth_km_;  % Updated to match the table column

%Read UMD info from Hua data
HuaMLD = readtable(filename3, 'Delimiter', ' ');
latitude_Hua_MLD = HuaMLD.latitude;
longitude_Hua_MLD = HuaMLD.longitude;
depth_Hua_MLD = HuaMLD.NVG_depth_km_;  

%Read UMD info from Hua data
HuaPVG = readtable(filename4, 'Delimiter', '\t');
latitude_Hua_PVG = HuaPVG.latitude;
longitude_Hua_PVG = HuaPVG.longitude;
depth_Hua_PVG = HuaPVG.PVG_150_depth_km;

%Read UMD info from Hopper
HopperLAB = xlsread(filename1, 2);
latitude_Hopper_LAB = HopperLAB(:,1);
longitude_Hopper_LAB = HopperLAB(:,2);
depth_Hopper_LAB = HopperLAB(:,3);

%Read UMD infor from Liu
LiuShallowMLD = readtable(filename5, 'Delimiter', '\t');
latitude_Liu_Shallow_MLD = LiuShallowMLD.Latitude;
longitude_Liu_Shallow_MLD = LiuShallowMLD.Longitude;
depth_Liu_Shallow_MLD = LiuShallowMLD.Depth;

LiuLAB = readtable(filename6, 'Delimiter', '\t');
latitude_Liu_LAB = LiuLAB{:,1};
longitude_Liu_LAB = LiuLAB{:,2};
depth_Liu_LAB = LiuLAB{:,3};

LiuDeepMLD = readtable('Liu_deep_MLD.txt', 'Delimiter', ' ');
latitude_Liu_Deep_MLD = LiuDeepMLD.Latitude;
longitude_Liu_Deep_MLD = LiuDeepMLD.Longitude;
depth_Liu_Deep_MLD = LiuDeepMLD.Depth;

% Read UMD info from Kreuger
Kreuger = readtable('Kreuger_UMD.txt', 'Delimiter', ',');
kreuger_lat = Kreuger.Lat;
kreuger_long = Kreuger.Long;
kreuger_depth = Kreuger.Depth;
kreuger_interp = Kreuger.Interp;

% Load Velocity model data from Schultz and Pelkum
Schultz_Pelkum_data = readtable('Schultz_Pelkum_vpvs.csv');
schultz_lat = Schultz_Pelkum_data.latitude;
schultz_long = Schultz_Pelkum_data.longitude;
schultz_depth = Schultz_Pelkum_data.depth;
schultz_vsv = Schultz_Pelkum_data.vsv;
schultz_vp = Schultz_Pelkum_data.vp;

%Get Moho data from Schmmandt
fid = fopen('US_Crust_GRL_2015_CrustThickness.txt', 'r');
Moho_data = textscan(fid, '%f%f%f', 'Delimiter', {' ', '\n'}, 'MultipleDelimsAsOne', true);
% Close the file
fclose(fid);
moho_lat = Moho_data{1};
moho_long = Moho_data{2};
moho_depth = Moho_data{3};

%% Split data into UMD1, UMD2 and UMD3
createTable = @(lat, long, depth, study, umd) [array2table([lat, long, depth], ...
    'VariableNames', {'Latitude', 'Longitude', 'Depth'}), ...
    table(repmat({study}, length(lat), 1), repmat({umd}, length(lat), 1), ...
    'VariableNames', {'Study', 'UMD'})];

% UMD1: Depth < 110 AND not interpreted as PVG
% Kreuger data with depth < 110 and not PVG
kreuger_UMD1 = ~(strcmp(kreuger_interp, 'PVG')) & kreuger_depth < 110;
UMD1_Hopper = createTable(latitude_Hopper_MLD(depth_Hopper_MLD < 110), longitude_Hopper_MLD(depth_Hopper_MLD < 110), depth_Hopper_MLD(depth_Hopper_MLD < 110), 'Hopper et al. 2018', 'UMD1');
UMD1_Abt = createTable(latitude_Abt(depth_Abt < 110), longitude_Abt(depth_Abt < 110), depth_Abt(depth_Abt < 110), 'Abt et al. 2010', 'UMD1');
UMD1_Hua = createTable(latitude_Hua_MLD(depth_Hua_MLD < 110), longitude_Hua_MLD(depth_Hua_MLD < 110), depth_Hua_MLD(depth_Hua_MLD < 110), 'Hua et al. 2023', 'UMD1');
UMD1_Liu_Shallow = createTable(latitude_Liu_Shallow_MLD(depth_Liu_Shallow_MLD < 110), longitude_Liu_Shallow_MLD(depth_Liu_Shallow_MLD < 110), depth_Liu_Shallow_MLD(depth_Liu_Shallow_MLD < 110), 'Liu and Shearer 2021', 'UMD1');
UMD1_Liu_LAB = createTable(latitude_Liu_LAB(depth_Liu_LAB < 110), longitude_Liu_LAB(depth_Liu_LAB < 110), depth_Liu_LAB(depth_Liu_LAB < 110), 'Liu and Shearer 2021', 'UMD1'); % Newly added Liu_LAB data for UMD1
UMD1_Kreuger = createTable(kreuger_lat(kreuger_UMD1), kreuger_long(kreuger_UMD1), kreuger_depth(kreuger_UMD1), 'Kreuger et al. 2021', 'UMD1');

% UMD2: Interp is PVG
UMD2_Hua = createTable(latitude_Hua_PVG, longitude_Hua_PVG, depth_Hua_PVG, 'Hua et al. 2023', 'UMD2');
UMD2_Kreuger = createTable(kreuger_lat(strcmp(kreuger_interp, 'PVG')), kreuger_long(strcmp(kreuger_interp, 'PVG')), kreuger_depth(strcmp(kreuger_interp, 'PVG')), 'Kreuger et al. 2021', 'UMD2');

% UMD3: Depth > 110 AND not interpreted as PVG
% Kreuger data with depth > 110 and not PVG
kreuger_UMD3 = ~(strcmp(kreuger_interp, 'PVG')) & kreuger_depth > 110;

UMD3_Hopper = createTable(latitude_Hopper_LAB(depth_Hopper_LAB > 110), longitude_Hopper_LAB(depth_Hopper_LAB > 110), depth_Hopper_LAB(depth_Hopper_LAB > 110), 'Hopper et al. 2018', 'UMD3');
UMD3_Abt = createTable(latitude_Abt(depth_Abt > 110), longitude_Abt(depth_Abt > 110), depth_Abt(depth_Abt > 110), 'Abt et al. 2010', 'UMD3');
UMD3_Hua = createTable(latitude_Hua_MLD(depth_Hua_MLD > 110), longitude_Hua_MLD(depth_Hua_MLD > 110), depth_Hua_MLD(depth_Hua_MLD > 110), 'Hua et al. 2023', 'UMD3');
UMD3_Liu = createTable(latitude_Liu_Deep_MLD(depth_Liu_Deep_MLD > 110), longitude_Liu_Deep_MLD(depth_Liu_Deep_MLD > 110), depth_Liu_Deep_MLD(depth_Liu_Deep_MLD > 110), 'Liu and Shearer 2021', 'UMD3');
UMD3_Kreuger = createTable(kreuger_lat(kreuger_UMD3), kreuger_long(kreuger_UMD3), kreuger_depth(kreuger_UMD3), 'Kreuger et al. 2021', 'UMD3');

UMD1_data = [UMD1_Hopper; UMD1_Abt; UMD1_Hua; UMD1_Liu_Shallow; UMD1_Liu_LAB; UMD1_Kreuger]; 
UMD2_data = [UMD2_Hua; UMD2_Kreuger];
UMD3_data = [UMD3_Hopper; UMD3_Abt; UMD3_Hua; UMD3_Liu; UMD3_Kreuger];
%%
rayP = linspace(0.02, 0.08, 60);

% Loop through each Moho depth and calculate average arrival times for Moho multiples (tPs, tPps, and tPss) 
for i = 1:length(moho_depth)
    % Find the closest point in the Schultz and Pelkum dataset based on geographic location
    distances = sqrt((moho_lat(i) - schultz_lat).^2 + (moho_long(i) - schultz_long).^2);
    [~, closest_idx] = min(distances);
    
    % Use the Vp and Vs values from the closest Schultz and Pelkum location
    Vp = schultz_vp(closest_idx);
    Vs = schultz_vsv(closest_idx);
    
    H = [moho_depth(i), 0];
    [tPs_all, tPps_all, tPss_all] = travelTimesAppx(Vp, Vs, H, rayP, 1, 2);
    
    % Averaging the arrival times across all rayP values
    tPs_avg(i) = mean(tPs_all);
    tPps_avg(i) = mean(tPps_all);
    tPss_avg(i) = mean(tPss_all);
end

% Initialize arrays to hold converted times
converted_times_UMD1 = zeros(size(UMD1_data.Depth));
converted_times_UMD2 = zeros(size(UMD2_data.Depth));
converted_times_UMD3 = zeros(size(UMD3_data.Depth));

% Convert UMD1 depths to arrival times
for i = 1:length(UMD1_data.Depth)
    [Vp, Vs] = find_closest_vp_vs(UMD1_data.Latitude(i), UMD1_data.Longitude(i), Schultz_Pelkum_data);
    convert_depth_to_time = @(depth) depth .* ((1/Vs) - (1/Vp));
    converted_times_UMD1(i) = convert_depth_to_time(UMD1_data.Depth(i));
end

% Convert UMD2 depths to arrival times
for i = 1:length(UMD2_data.Depth)
    [Vp, Vs] = find_closest_vp_vs(UMD2_data.Latitude(i), UMD2_data.Longitude(i), Schultz_Pelkum_data);
    convert_depth_to_time = @(depth) depth .* ((1/Vs) - (1/Vp));
    converted_times_UMD2(i) = convert_depth_to_time(UMD2_data.Depth(i));
end

% Convert UMD3 depths to arrival times
for i = 1:length(UMD3_data.Depth)
    [Vp, Vs] = find_closest_vp_vs(UMD3_data.Latitude(i), UMD3_data.Longitude(i), Schultz_Pelkum_data);
    convert_depth_to_time = @(depth) depth .* ((1/Vs) - (1/Vp));
    converted_times_UMD3(i) = convert_depth_to_time(UMD3_data.Depth(i));
end

%%
% Find matching Moho depths
matching_moho_depths_UMD1 = find_matching_moho_depths(UMD1_data.Latitude, UMD1_data.Longitude, moho_lat, moho_long, moho_depth);
matching_moho_depths_UMD2 = find_matching_moho_depths(UMD2_data.Latitude, UMD2_data.Longitude, moho_lat, moho_long, moho_depth);
matching_moho_depths_UMD3 = find_matching_moho_depths(UMD3_data.Latitude, UMD3_data.Longitude, moho_lat, moho_long, moho_depth);

moho_depth = moho_depth(:);
tPps_avg = tPps_avg(:);
tPss_avg = tPss_avg(:);

% Create a grid over the range of your data
[x, y] = meshgrid(linspace(min(moho_depth), max(moho_depth), 100), ...
                  linspace(min(min(tPps_avg), min(tPss_avg)), max(max(tPps_avg), max(tPss_avg)), 100));

% Flatten the grid matrices for ksdensity input
x_flat = x(:);
y_flat_tPps = y(:);
y_flat_tPss = y(:);

% Perform two dimensional kernel density estimation for tPps_avg
[f_tPps, ~] = ksdensity([moho_depth, tPps_avg], [x_flat, y_flat_tPps],  'Function', 'pdf');

% Perform kernel density estimation for tPss_avg
[f_tPss, ~] = ksdensity([moho_depth, tPss_avg], [x_flat, y_flat_tPss], 'Function', 'pdf');

% Reshape the output for contour plotting
z_tPps = reshape(f_tPps, size(x));
z_tPss = reshape(f_tPss, size(x));
%%
figure(1);clf;

% Plot density estimation for tPps_avg
contour(x, y, z_tPps, 50, 'LineColor', 'r');  % 50 contour levels, adjust as needed
hold on;
% Plot density estimation for tPss_avg
contour(x, y, z_tPss, 50, 'LineColor', 'r');  % 50 contour levels, adjust as needed

% Define the datasets and matching moho depths arrays
matching_moho_depths = {matching_moho_depths_UMD1, matching_moho_depths_UMD2, matching_moho_depths_UMD3};
converted_times = {converted_times_UMD1, converted_times_UMD2, converted_times_UMD3};

datasets = {UMD1_data, UMD2_data, UMD3_data}
% Define marker properties for each study and UMD
markerProperties(1) = struct('Study', 'Hopper et al. 2018', 'UMD', 'UMD1', 'markerSize', 20, 'markerAlpha', 0.5, 'edgeColor', 'b');
markerProperties(2) = struct('Study', 'Hopper et al. 2018', 'UMD', 'UMD3', 'markerSize', 20, 'markerAlpha', 0.5, 'edgeColor', 'b');

markerProperties(3) = struct('Study', 'Abt et al. 2010', 'UMD', 'UMD1', 'markerSize', 60, 'markerAlpha', 1, 'edgeColor', 'k');
markerProperties(4) = struct('Study', 'Abt et al. 2010', 'UMD', 'UMD3', 'markerSize', 60, 'markerAlpha', 1, 'edgeColor', 'k');

markerProperties(5) = struct('Study', 'Hua et al. 2023', 'UMD', 'UMD1', 'markerSize', 50, 'markerAlpha', 1, 'edgeColor', 'k');
markerProperties(6) = struct('Study', 'Hua et al. 2023', 'UMD', 'UMD2', 'markerSize', 70, 'markerAlpha', 1, 'edgeColor', 'k');
markerProperties(7) = struct('Study', 'Hua et al. 2023', 'UMD', 'UMD3', 'markerSize', 70, 'markerAlpha', 1, 'edgeColor', 'k');

markerProperties(8) = struct('Study', 'Kreuger et al. 2021', 'UMD', 'UMD1', 'markerSize', 50, 'markerAlpha', 1, 'edgeColor', 'k');
markerProperties(9) = struct('Study', 'Kreuger et al. 2021', 'UMD', 'UMD2', 'markerSize', 70, 'markerAlpha', 1, 'edgeColor', 'k');
markerProperties(10) = struct('Study', 'Kreuger et al. 2021', 'UMD', 'UMD3', 'markerSize', 70, 'markerAlpha', 1, 'edgeColor', 'k');

markerProperties(11) = struct('Study', 'Liu and Shearer 2021', 'UMD', 'UMD1', 'markerSize', 50, 'markerAlpha', 1, 'edgeColor', 'k');
markerProperties(12) = struct('Study', 'Liu and Shearer 2021', 'UMD', 'UMD3', 'markerSize', 70, 'markerAlpha', 1, 'edgeColor', 'k');

tic

% Initialize containers for legend handles and labels
legendHandles = [];
legendLabels = {};

hold on;  % Move hold on outside the loop

for d = 1:length(datasets)
    data = datasets{d};
    disp(['Processing dataset ', num2str(d)]);  % Debugging output
    % Create a unique identifier for each group of data points
    groupID = strcat(cellstr(data.Study), '_', cellstr(data.UMD));

    uniqueGroups = unique(groupID);
    
    for ug = 1:length(uniqueGroups)
        idx = strcmp(groupID, uniqueGroups{ug});
        disp(['Processing unique group ', uniqueGroups{ug}]);  % Debugging output
        % Parse the unique group identifier to get the study and UMD
        tokens = split(uniqueGroups{ug}, '_');
        study = tokens{1};
        umd = tokens{2};
        
        % Look up marker properties for the current study and UMD
        properties = [];
        for mp = 1:length(markerProperties)
            if strcmp(markerProperties(mp).Study, study) && strcmp(markerProperties(mp).UMD, umd)
                properties = markerProperties(mp);
                break;
            end
        end
        
        if ~isempty(properties)
            markerSize = properties.markerSize;
            markerAlpha = properties.markerAlpha;
            edgeColor = properties.edgeColor;
            
            marker = get_marker(study);
            color = get_color(umd);
            
            % Use scatter function to plot with specified properties
            scatter(matching_moho_depths{d}(idx), converted_times{d}(idx), markerSize, ...
                'Marker', marker, 'MarkerFaceColor', color, ...
                'MarkerEdgeColor', edgeColor, 'MarkerFaceAlpha', markerAlpha);
        end
    end
end

%for Kind et al. 2020 and Liu and Shearer 2021 data
% Convert these depth ranges to time ranges
kind_depths = [150 200];  % LABc

convert_depth_to_time = @(depth) depth .* ((1/4.3) - (1/8));  % Assuming the same Vp and Vs values

kind_times = convert_depth_to_time(kind_depths);


% Assuming the average Moho depth for Kind et al. data is 35 km
moho_depth_kind = 40;
errorbar(moho_depth_kind, mean(kind_times(1, :)), diff(kind_times(1, :))/2, 'h', 'MarkerSize', 15, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'Color', 'k','LineWidth',1.5);

% After your loop, create legend entries for each study
uniqueStudies = unique({markerProperties.Study});
legendSymbolSize = 100;  % Adjust this value to your desired size

for i = 1:length(uniqueStudies)
    % Find the properties for this study (assuming each study uses the same marker)
    for mp = 1:length(markerProperties)
        if strcmp(markerProperties(mp).Study, uniqueStudies{i})
            properties = markerProperties(mp);
            break;
        end
    end
    legendHandles(end+1) = scatter(NaN, NaN, legendSymbolSize, ...
        'Marker', get_marker(properties.Study), ...
        'MarkerEdgeColor', properties.edgeColor);
    legendLabels{end+1} = [uniqueStudies{i}, ' ', get_publication_year(uniqueStudies{i})]; %#ok<SAGROW>
end

% Create the legend
legendHandles(end+1) = scatter(NaN, NaN, legendSymbolSize, 'h', 'MarkerEdgeColor', 'k');  % Kind et al. 2020
legendLabels{end+1} = 'Kind et al. 2020';

legendHandles(end+1) = scatter(NaN, NaN, legendSymbolSize, 'v', 'MarkerEdgeColor', 'k');  % Liu and Shearer 2021 MLD
legendLabels{end+1} = 'Liu and Shearer 2021';

legend(legendHandles, legendLabels, 'Location', 'NorthWest','FontSize',16);

ylabel('Delay Time [s]','FontSize',35)
xlabel('Moho Depth Z_c [km]','FontSize',35) 
ax = gca;
ax.XAxis.FontSize = 30;  
ax.YAxis.FontSize = 30;
%print(figure(1),'./Figures/Figure_1/figure1a','-vector','-dpdf','-r0');
%%
% Extract the depth data from each dataset
depths_UMD1 = UMD1_data.Depth;  % Assuming Depth is the name of the depth column
depths_UMD2 = UMD2_data.Depth;
depths_UMD3 = UMD3_data.Depth;

% Histogram settings
binWidth = 5;  % Modify this according to the precision you want

% Plot for UMD1
figure(2);clf; % Open a new figure for UMD1
histogram(depths_UMD1, 'BinWidth', binWidth, 'FaceColor', 'blue', 'EdgeColor', 'blue');
title('Histogram of Depths for UMD1');
xlabel('Depth [km]');
ylabel('Frequency');
legend('UMD 1')
%print(figure(2),'./Figures/Figure_1/figure1b','-vector','-dpdf','-r0');

% Plot for UMD2
figure(3);clf; % Open a new figure for UMD2
histogram(depths_UMD2, 'BinWidth', binWidth, 'FaceColor', 'magenta', 'EdgeColor', 'magenta');
title('Histogram of Depths for UMD2');
xlabel('Depth [km]');
ylabel('Frequency');
legend('UMD 2')
%print(figure(3),'./Figures/Figure_1/figure1c','-vector','-dpdf','-r0');

% Plot for UMD3
figure(4);clf % Open a new figure for UMD3
histogram(depths_UMD3, 'BinWidth', binWidth, 'FaceColor', 'green', 'EdgeColor', 'green');
title('Histogram of Depths for UMD3');
xlabel('Depth [km]');
ylabel('Frequency');
legend('UMD 3')
%print(figure(4),'./Figures/Figure_1/figure1d','-vector','-dpdf','-r0');
%%
% Initialize the contaminated data structure
contaminated = struct('Longitude', [], 'Latitude', [], 'Study', [], 'Depth', []);

% Filter for UMD1
idx_UMD1 = (converted_times_UMD1 >= 8 & converted_times_UMD1 <= 12) & ...
           (matching_moho_depths_UMD1 >= 20 & matching_moho_depths_UMD1 <= 35);
contaminated.Longitude = [contaminated.Longitude; UMD1_data.Longitude(idx_UMD1)];
contaminated.Latitude = [contaminated.Latitude; UMD1_data.Latitude(idx_UMD1)];
contaminated.Study = [contaminated.Study; UMD1_data.Study(idx_UMD1)];
contaminated.Depth = [contaminated.Depth; UMD1_data.Depth(idx_UMD1)];

% Filter for UMD2
idx_UMD2 = (converted_times_UMD2 >= 10 & converted_times_UMD2 <= 17) & ...
           (matching_moho_depths_UMD2 >= 20 & matching_moho_depths_UMD2 <= 50);
contaminated.Longitude = [contaminated.Longitude; UMD2_data.Longitude(idx_UMD2)];
contaminated.Latitude = [contaminated.Latitude; UMD2_data.Latitude(idx_UMD2)];
contaminated.Study = [contaminated.Study; UMD2_data.Study(idx_UMD2)];
contaminated.Depth = [contaminated.Depth; UMD2_data.Depth(idx_UMD2)];

% Filter for UMD3
idx_UMD3 = (converted_times_UMD3 >= 10 & converted_times_UMD3 <= 22) & ...
           (matching_moho_depths_UMD3 >= 20 & matching_moho_depths_UMD3 <= 47.5);
contaminated.Longitude = [contaminated.Longitude; UMD3_data.Longitude(idx_UMD3)];
contaminated.Latitude = [contaminated.Latitude; UMD3_data.Latitude(idx_UMD3)];
contaminated.Study = [contaminated.Study; UMD3_data.Study(idx_UMD3)];
contaminated.Depth = [contaminated.Depth; UMD3_data.Depth(idx_UMD3)];
%%
% List of studies to loop over
studies = {'Hopper et al. 2018', 'Abt et al. 2010', 'Hua et al. 2023', 'Kreuger et al. 2021', 'Liu and Shearer 2021'};

% Initialize the map of the US with a Mercator projection
figure(5); clf;

% Create a map axes with slightly adjusted limits
ax = axesm('MapProjection', 'mercator', ...
           'MapLatLimit', [24 50], ...
           'MapLonLimit', [-126 -66]);

% Display a graticule on the map
framem on;

% Add latitude (parallel) labels on the left side of the map
plabel('PLabelLocation', 10, 'PLabelRound', -1, 'PLabelMeridian', 'west');

% Add longitude (meridian) labels on the bottom of the map
mlabel('MLabelLocation', 10, 'MLabelRound', -1, 'MLabelParallel', 'south');
hold on;

greyColor = [0.5 0.5 0.5]; % Defining grey color

% Latitude and Longitude boundaries for contiguous US
lat_bounds = [24.396308, 49.384358];
lon_bounds = [-125.001650, -66.934570];

for study = studies
    study = study{1};  % Extract string from cell
    marker = get_marker(study);
    markerSize = get_marker_size(study);
    
    idx = ~idx_UMD1 & strcmp(UMD1_data.Study, study) ...
        & UMD1_data.Latitude >= lat_bounds(1) & UMD1_data.Latitude <= lat_bounds(2) ...
        & UMD1_data.Longitude >= lon_bounds(1) & UMD1_data.Longitude <= lon_bounds(2);
    [x, y] = mfwdtran(UMD1_data.Latitude(idx), UMD1_data.Longitude(idx));
    scatter(x, y, markerSize, greyColor, 'o', 'filled', 'MarkerFaceAlpha', 0.2);
    
    idx = ~idx_UMD2 & strcmp(UMD2_data.Study, study) ...
        & UMD2_data.Latitude >= lat_bounds(1) & UMD2_data.Latitude <= lat_bounds(2) ...
        & UMD2_data.Longitude >= lon_bounds(1) & UMD2_data.Longitude <= lon_bounds(2);
    [x, y] = mfwdtran(UMD2_data.Latitude(idx), UMD2_data.Longitude(idx));
    scatter(x, y, markerSize, greyColor, 'o', 'filled', 'MarkerFaceAlpha', 0.2);
    
    idx = ~idx_UMD3 & strcmp(UMD3_data.Study, study) ...
        & UMD3_data.Latitude >= lat_bounds(1) & UMD3_data.Latitude <= lat_bounds(2) ...
        & UMD3_data.Longitude >= lon_bounds(1) & UMD3_data.Longitude <= lon_bounds(2);
    [x, y] = mfwdtran(UMD3_data.Latitude(idx), UMD3_data.Longitude(idx));
    scatter(x, y, markerSize, greyColor, 'o', 'filled', 'MarkerFaceAlpha', 0.2);
end

% Now plotting contaminated data points using geoshow to ensure they are on top
for study = studies
    study = study{1};  % Extract string from cell
    marker = get_marker(study);
    markerSize = get_marker_size(study);
    
    idx = strcmp(contaminated.Study, study);
    geoshow(contaminated.Latitude(idx), contaminated.Longitude(idx), 'DisplayType', 'point', 'Marker', marker, 'MarkerSize', markerSize, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r');
end

% Draw the physio boundaries after plotting all points
hold on;
geoshow('physio.shp', 'DisplayType', 'polygon', 'EdgeColor', 'k', 'LineWidth', 1, 'FaceColor', 'none');

blueColor = [0 0 1];
redColor = [1 0 0];

% Now, adding the legend for red and blue points
lgd = legend([scatter([],[],'o','filled', 'MarkerFaceColor', redColor), scatter([],[],'o','filled', 'MarkerFaceColor', blueColor)], ...
    {'Likely to suffer contamination by moho multiples', 'Unlikely to suffer contamination by moho multiples'});
set(lgd,'Location','southwest');
%print(figure(5),'./Figures/Figure_1/figure1e','-vector','-dpdf','-r0');


