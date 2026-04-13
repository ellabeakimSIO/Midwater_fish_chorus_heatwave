% -------------------------------------------------------------------------
% 01_plot_figure1a_study_region_map.m
%
% Purpose:
% Generate the main study-region map used for Figure 1A, showing UF440 chorus
% presence and median sunset SPL across sites in the Southern California Bight.
%
% In the paper:
% This script reproduces the main map component of Figure 1A in Kim et al.
% (2026). Circle size scales with the proportion of total recording hours
% containing chorus detections, and circle color represents median sunset
% chorus SPL. Additional symbols indicate sites with low or no observed UF440
% chorus presence.
%
% Inputs:
%   - GEBCO bathymetry netCDF
%   - CINMS sanctuary boundary shapefile
%   - Sites.csv
%   - UF440 site summary table(s) with:
%       Site, Presence_percent, Median_UF440_sunset_SPL
%
% Outputs:
%   - Figure1A_study_region_map.png
%   - Figure1A_study_region_map.svg
%
% Notes:
% - This script reproduces the Figure 1A map layer only.
% - Additional layout elements in the final assembled figure
%   (for example labels, inset box, or panel B assembly) were made in Adobe Illustrator
% -------------------------------------------------------------------------

clear; clc;

% =========================
% User settings
% =========================
dataDir = 'data';
outputDir = 'figures';

if ~exist(outputDir, 'dir')
    mkdir(outputDir);
end

ncfile = fullfile(dataDir, 'gebco_sc_bight.nc');
shpfile = fullfile(dataDir, 'cinms_boundary.shp');
sites_csv = fullfile(dataDir, 'Sites.csv');

summary_main_mat = fullfile(dataDir, 'UF440_site_summary_no20iso_noday.mat');
summary_main_csv = fullfile(dataDir, 'UF440_site_summary_no20iso_noday.csv');
summary_ci_mat   = fullfile(dataDir, 'UF440_site_summary_no20iso_noday_CI010305.mat');
summary_ci_csv   = fullfile(dataDir, 'UF440_site_summary_no20iso_noday_CI010305.csv');

outBase = fullfile(outputDir, 'Figure1A_study_region_map');

% =========================
% Load bathymetry
% =========================
lat = double(ncread(ncfile, 'lat'));
lon = double(ncread(ncfile, 'lon'));
[mlon, mlat] = meshgrid(lon, lat);

bathymetry = double(ncread(ncfile, 'elevation'));
bathymetry(bathymetry > 0) = 1000;

% =========================
% Sanctuary outline
% =========================
scinms = shaperead(shpfile, 'UseGeoCoords', true);

% =========================
% Read site coordinates
% =========================
Tsites = readtable(sites_csv, 'TextType', 'string');

% Keep only non-SDT sites plus SDT_HP
isSDT = startsWith(upper(strtrim(Tsites.Site)), "SDT_");
keepSDT_HP = startsWith(upper(strtrim(Tsites.Site)), "SDT_HP");
Tsites = Tsites(~isSDT | keepSDT_HP, :);

Tsites.lat_dd = arrayfun(@parseDegMin, Tsites.Lat);
Tsites.lon_dd = arrayfun(@parseDegMin, Tsites.Long);

% =========================
% UF440 bubble sites
% =========================
bubbleSites = ["CI01","CI03","CI05","South_P","South_T","North_B"];

bub_lat = nan(numel(bubbleSites),1);
bub_lon = nan(numel(bubbleSites),1);

for i = 1:numel(bubbleSites)
    s = bubbleSites(i);

    if s == "North_B"
        idx = startsWith(upper(strtrim(Tsites.Site)), "NORTH_B") | ...
              startsWith(upper(strtrim(Tsites.Site)), "CINMS_B");
    else
        idx = upper(strtrim(Tsites.Site)) == upper(s);
    end

    bub_lat(i) = mean(Tsites.lat_dd(idx), 'omitnan');
    bub_lon(i) = mean(Tsites.lon_dd(idx), 'omitnan');
end

% SDT_HP coordinates (red triangle)
idxHP = startsWith(upper(strtrim(Tsites.Site)), "SDT_HP");
sdt_hp_lat = mean(Tsites.lat_dd(idxHP), 'omitnan');
sdt_hp_lon = mean(Tsites.lon_dd(idxHP), 'omitnan');

% =========================
% Load UF440 summary stats
% =========================
summaryAll = loadSummary(summary_main_mat, summary_main_csv);
summaryCI  = loadSummary(summary_ci_mat, summary_ci_csv);
summaryTbl = [summaryAll; summaryCI];

% Normalize older names if needed
summaryTbl.Site = replace(summaryTbl.Site, "CINMS_B", "North_B");
summaryTbl.Site = replace(summaryTbl.Site, "SOCAL_P", "South_P");
summaryTbl.Site = replace(summaryTbl.Site, "SOCAL_T", "South_T");

[~, ia] = unique(string(summaryTbl.Site), 'stable');
summaryTbl = summaryTbl(ia, :);

Presence_percent = nan(numel(bubbleSites),1);
Median_SPL = nan(numel(bubbleSites),1);

for i = 1:numel(bubbleSites)
    s = bubbleSites(i);
    idx = string(summaryTbl.Site) == s;
    if any(idx)
        j = find(idx, 1, 'first');
        Presence_percent(i) = summaryTbl.Presence_percent(j);
        Median_SPL(i) = summaryTbl.Median_UF440_sunset_SPL(j);
    end
end

Presence_round = round(Presence_percent);

% =========================
% Bathymetry styling
% =========================
bathyCol = bone(15);
bathyCol = bathyCol(1:13, :);
bathyCol(14, :) = [1 1 1];

bathy_min = min(bathymetry(bathymetry <= 0));
lvls = linspace(bathy_min, 0, 14);

% =========================
% Figure sizing
% =========================
W_in = 7.5;
H_in = 3.7;

fig = figure(3); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 W_in H_in]);
set(fig, 'PaperUnits', 'inches', 'PaperSize', [W_in H_in], ...
         'PaperPosition', [0 0 W_in H_in], 'PaperPositionMode', 'manual');

axPos = [0.08 0.12 0.79 0.82];

% =========================
% Axis 1: bathymetry + sanctuary outline
% =========================
ax1 = axes(fig, 'Position', axPos);
[~, hc] = contourf(ax1, mlon, mlat, bathymetry.', lvls);
set(ax1, 'YDir', 'normal');
hc.LineColor = [0.2 0.2 0.2];
colormap(ax1, bathyCol);
caxis(ax1, [lvls(1) lvls(end)]);
hold(ax1, 'on');

geoshow(ax1, scinms, 'FaceColor', 'none', 'EdgeColor', 'm', 'LineWidth', 2);

xlim(ax1, [min(lon) max(lon)]);
ylim(ax1, [min(lat) max(lat)]);
xlabel(ax1, '');
ylabel(ax1, '');
box(ax1, 'on');

% =========================
% Axis 2: symbols + bubbles
% =========================
ax2 = axes(fig, 'Position', axPos);
set(ax2, 'Color', 'none', ...
         'XLim', get(ax1, 'XLim'), ...
         'YLim', get(ax1, 'YLim'), ...
         'XTick', [], 'YTick', [], ...
         'Box', 'off');
hold(ax2, 'on');

% Triangle sites
isTriSites = ismember(lower(strtrim(Tsites.Color)), ["black","gray","grey","white"]);

tri_lon = Tsites.lon_dd(isTriSites);
tri_lat = Tsites.lat_dd(isTriSites);
tri_clr = lower(strtrim(string(Tsites.Color(isTriSites))));

tri_rgb = zeros(sum(isTriSites), 3);
for ii = 1:numel(tri_clr)
    tri_rgb(ii,:) = colorToRGB(tri_clr(ii));
end

triSize = 55;
scatter(ax2, tri_lon, tri_lat, triSize, tri_rgb, '^', ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.8);

% SDT_HP as red triangle
hpSize = 65;
scatter(ax2, sdt_hp_lon, sdt_hp_lat, hpSize, [1 0 0], '^', ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 0.9);

% =========================
% Bubble sizing
% =========================
minA = 60;
maxA = 700;

pp = Presence_round;
pp(isnan(pp)) = 0;
pp(pp < 0) = 0;

pp_min = min(pp);
pp_max = max(pp);

if pp_max == pp_min
    bubbleA = repmat((minA + maxA) / 2, size(pp));
else
    bubbleA = minA + (sqrt(pp) - sqrt(pp_min)) ./ (sqrt(pp_max) - sqrt(pp_min)) .* (maxA - minA);
end

% =========================
% Bubble colors: discrete SPL bins
% =========================
spl = Median_SPL;

nbins = 6;
bin_edges = 64:2:76;
bin_labs = 64:2:74;

spl_clip = spl;
spl_clip = max(spl_clip, bin_edges(1));
spl_clip = min(spl_clip, bin_edges(end) - eps);

bin_idx = discretize(spl_clip, bin_edges);

scatter(ax2, bub_lon, bub_lat, bubbleA, bin_idx, ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);

cmap7 = jet(7);
cmap6 = cmap7([1 3:7], :);
colormap(ax2, cmap6);
caxis(ax2, [0.5 nbins + 0.5]);

% SPL colorbar
cb = colorbar(ax2, 'eastoutside');
cb.TickDirection = 'out';
cb.Label.String = 'Median UF440 sunset SPL (dB re 1 \muPa)';
cb.Ticks = 1:nbins;
cb.TickLabels = string(bin_labs);

cbW = 0.022;
gap = 0.008;
cbX = axPos(1) + axPos(3) + gap;
cb.Position = [cbX 0.20 cbW 0.60];

% =========================
% Save outputs
% =========================
pngFile = strcat(outBase, '.png');
svgFile = strcat(outBase, '.svg');

exportgraphics(fig, pngFile, 'Resolution', 300);
set(fig, 'Renderer', 'painters');
print(fig, svgFile, '-dsvg', '-painters');

fprintf('Saved:\n  %s\n  %s\n', pngFile, svgFile);

% =========================
% Helper functions
% =========================
function dd = parseDegMin(s)
    s = string(s); s = strtrim(s);
    parts = split(s);
    hemi = upper(parts(end));
    dm = parts(1);
    dmParts = split(dm, "-");
    deg = str2double(dmParts(1));
    minutes = str2double(dmParts(2));
    dd = deg + minutes / 60;
    if hemi == "S" || hemi == "W"
        dd = -dd;
    end
end

function summaryTbl = loadSummary(matFile, csvFile)
    if exist(matFile, 'file') == 2
        S = load(matFile);
        if isfield(S, 'summaryTbl')
            summaryTbl = S.summaryTbl;
        else
            error('MAT file %s does not contain variable summaryTbl.', matFile);
        end
    elseif exist(csvFile, 'file') == 2
        summaryTbl = readtable(csvFile, 'TextType', 'string');
    else
        error('Could not find %s or %s', matFile, csvFile);
    end

    req = ["Site","Presence_percent","Median_UF440_sunset_SPL"];
    for r = req
        if ~ismember(r, string(summaryTbl.Properties.VariableNames))
            error('Summary file missing required column: %s', r);
        end
    end

    summaryTbl.Site = string(summaryTbl.Site);
end

function rgb = colorToRGB(c)
    c = lower(strtrim(string(c)));
    switch c
        case "red"
            rgb = [1 0 0];
        case "orange"
            rgb = [1 0.5 0];
        case {"gray","grey"}
            rgb = [0.5 0.5 0.5];
        case "black"
            rgb = [0 0 0];
        case "white"
            rgb = [1 1 1];
        otherwise
            rgb = [0 0 1];
    end
end
