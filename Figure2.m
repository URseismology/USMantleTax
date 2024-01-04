close all; clear;clc;
localBaseDir = '/scratch/tolugboj_lab/';
%localBaseDir = '/Users/stevecarr/Documents/bluehive';
pathname= [localBaseDir '/Prj7_RadonT/Prj7_US_Earthscope/2_Data/'];
addpath([localBaseDir '/Prj7_RadonT/Prj7_US_Earthscope/figures_4pub/m_map/']);
addpath([localBaseDir '/Prj7_RadonT/Prj7_US_Earthscope/figures_4pub/'])
FUNCDIR = [localBaseDir '/Prj7_RadonT/RFImager_EvansVersion/1_Functions/'];

addpath(FUNCDIR); addpath(pathname);
RFDIR = [localBaseDir 'Prj7_RadonT/2_Data/Parallel_Download/MTCRF/'];
%%
% Map Projection
hh = figure (1);
clf;
fs=10;
m_proj('Mercator','lon',[-128 -63],'lat',[25 49.5]); %continental US
shading flat;

%Geologic Ages data
age_data = 'global-ages-0705-1x1.xyz.txt';
data = readtable(age_data, 'Delimiter', '\t', 'Format', '%f%f%f');
lon = data{:, 1};
lat = data{:, 2};
age = data{:, 3};
[lon2, lat2, age2] = xyz2grid(lon, lat, age);
m_contourf(lon2, lat2, age2, 'linecolor', 'none');

cdata = [...
    0 245 236 191 540 245 236 191;...
    540 240 184 116 2500 240 184 116;...
    2500 241 202 183 3500 241 202 183];
dlmwrite('mycmap.cpt', cdata, ' ');
cptcmap('mycmap', 'mapping', 'direct', 'ncol', 30*5);

%%

% Read the data from the acceptedStations_nocrust2.txt file
filename = 'acceptedStations_nocrust2.txt';
stationsData = readtable(filename);

% Separate the data into TA and non-TA stations
TA_stations = stationsData(strcmp(stationsData.Network, 'TA'), :);
non_TA_stations = stationsData(~strcmp(stationsData.Network, 'TA'), :);

% Extract longitude and latitude for TA and non-TA stations
lon_TA = TA_stations.faceLon;
lat_TA = TA_stations.faceLat;
lon_non_TA = non_TA_stations.faceLon;
lat_non_TA = non_TA_stations.faceLat;

% Plot the TA stations as blue
h_TA = m_line(lon_TA, lat_TA, 'linestyle', 'none', 'marker', '^', 'color', 'k', 'MarkerFaceColor', 'b','MarkerSize',8);

% Plot the non-TA stations as red
h_non_TA = m_line(lon_non_TA, lat_non_TA, 'linestyle', 'none', 'marker', '^', 'color', 'k', 'MarkerFaceColor', 'g','MarkerSize',8);

% Highlight two specific stations
% Coordinates for the two specific stations
lon_highlight = [-113.9; -68.2];
lat_highlight = [46.8; 44.7];

% Plot the highlighted stations as red triangles and make them slightly larger
h_highlight = m_line(lon_highlight, lat_highlight, 'linestyle', 'none', 'marker', '^', 'color', 'k', 'MarkerFaceColor', 'r', 'MarkerSize', 10);

% Add title and potentially other annotations
title('Seismic Stations Locations');

%%
% Read the shapefile using MATLAB's shaperead
S = shaperead('physio.shp');

% Loop through the shapes and plot them on the m_map
for k = 1:length(S)
    m_line(S(k).X, S(k).Y, 'color', 'k', 'LineWidth', 1);
end

% Coastlines and Grids
m_coast('line','color','k','linewidth',1);
m_grid('tickdir','out','linewidth',3,'fontsize',16);

% Add Legend
lgd = legend([h_TA(1), h_non_TA(1), h_highlight(1)], {'USArray TA stations', 'Other stations', 'Example stations'}, 'Location', 'southeast');
colorbar('southoutside')