close all; clear;clc;
addpath(FUNCDIR); addpath(pathname);
RFDIR = [localBaseDir 'Prj7_RadonT/2_Data/Parallel_Download/MTCRF/'];
%%
% Map Projection
hh = figure (1);
clf;
fs=10;
m_proj('Mercator','lon',[-126 -63],'lat',[25 49.5]); %continental US
shading flat;

[CS,CH]=m_etopo2('contourf',[-5000:500:0 250:250:3000],'edgecolor','none');
colormap([ m_colmap('blues',80); flipud( gray(48) ) ]);
colorbar('southoutside')

%%

% Read the data from the acceptedStations_nocrust2.txt file
filename = 'acceptedStations_nocrust2.txt';
stationsData = readtable(filename);

% Separate the data into TA and non-TA stations
TA_stations = stationsData(strcmp(stationsData.Network, 'TA'), :);
non_TA_stations = stationsData(~strcmp(stationsData.Network, 'TA'), :);

% Extract longitude and latitude for TA and non-TA stations
lon_TA = TA_stations.Lon;
lat_TA = TA_stations.Lat;
lon_non_TA = non_TA_stations.Lon;
lat_non_TA = non_TA_stations.Lat;

% Plot the TA stations as blue
h_TA = m_line(lon_TA, lat_TA, 'linestyle', 'none', 'marker', '^', 'color', 'k', 'MarkerFaceColor', 'g','MarkerSize',8);

% Plot the non-TA stations as red
h_non_TA = m_line(lon_non_TA, lat_non_TA, 'linestyle', 'none', 'marker', '^', 'color', 'k', 'MarkerFaceColor', 'b','MarkerSize',8);

% Highlight two specific stations
% Coordinates for the two specific stations
lon_highlight = [-113.9; -68.2];
lat_highlight = [46.8; 44.7];

% Plot the highlighted stations as red triangles and make them slightly larger
h_highlight = m_line(lon_highlight, lat_highlight, 'linestyle', 'none', 'marker', '^', 'color', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);

% Add title and potentially other annotations
title('Seismic Stations Locations');

%%
% Read US geological provinces
S = shaperead('physio.shp');

% List of provinces to plot
provincesToPlot = {'COASTAL PLAIN','PACIFIC BORDER'};

% Loop through the shapes and plot only the specified provinces
for k = 1:length(S)
    if any(strcmp(S(k).PROVINCE, provincesToPlot))
        m_line(S(k).X, S(k).Y, 'color', 'r', 'LineWidth', 3);
    end
end

% Read the boundary coordinates from the CSV file
Basin_and_range = readtable('Basin_and_range.csv');
Appalatians = readtable('Appalatians.csv');
Snake_river_plain = readtable('Snake_river_plain.csv');
Colorado_plateau = readtable('Colorado_plateau.csv');
Rocky_mountains = readtable('Rocky_mountains.csv');
Archean_boundary = readtable('Archean_boundary.csv');

% Overlay the boundary as a red dashed line
m_line(Basin_and_range.Longitude, Basin_and_range.Latitude, 'color', 'r', 'LineStyle', '-', 'LineWidth', 3);
m_line(Appalatians.Longitude, Appalatians.Latitude, 'color', 'r', 'LineStyle', '-', 'LineWidth', 3);
m_line(Snake_river_plain.Longitude, Snake_river_plain.Latitude, 'color', 'r', 'LineStyle', '-', 'LineWidth', 3);
m_line(Colorado_plateau.Longitude, Colorado_plateau.Latitude, 'color', 'r', 'LineStyle', '-', 'LineWidth', 3);
m_line(Rocky_mountains.Longitude, Rocky_mountains.Latitude, 'color', 'r', 'LineStyle', '-', 'LineWidth', 3);
m_line(Archean_boundary.Longitude, Archean_boundary.Latitude, 'color', 'm', 'LineStyle', '--', 'LineWidth', 6);

%%
% Coastlines and Grids
m_coast('line','color','k','linewidth',1);
m_grid('tickdir','out','linewidth',3,'fontsize',16);

% % states
M=m_shaperead('ne_50m_admin_1_states_provinces');
for k=1:length(M.ncst)
    m_line(M.ncst{k}(:,1),M.ncst{k}(:,2),'color','k','linewi',0.1);
end
