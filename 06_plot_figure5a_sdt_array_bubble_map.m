% -------------------------------------------------------------------------
% 06_plot_figure5a_sdt_array_bubble_map.m
%
% Purpose:
% Generate the San Diego regional array bubble map used for Figure 5A,
% showing UF440 chorus presence and median sunset SPL across SDT sites and
% South_T on a GEBCO bathymetry background.
%
% In the paper:
% This script reproduces the core site-bubble map component of Figure 5A in
% Kim et al. (2026), where circle size scales with the proportion of total
% recording hours containing chorus detections and circle color represents
% median sunset chorus SPL.
%
% Inputs:
%   - GEBCO bathymetry netCDF
%   - Sites.csv
%   - figure5a_site_summary.csv (or .mat), containing:
%       Site, Presence_percent, Median_UF440_sunset_SPL
%
% Outputs:
%   - Figure5A_SDT_array_bubble_map.png
%   - Figure5A_SDT_array_bubble_map.svg
%
% Notes:
% - This script reproduces the map layer for quantified SDT sites and South_T.
% - Additional labels and layout edits were added downstream in
%   Illustrator for the final assembled figure.
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

ncfile = fullfile(dataDir, 'gebco_sdt_region.nc');
sites_csv = fullfile(dataDir, 'Sites.csv');
summary_mat = fullfile(dataDir, 'figure5a_site_summary.mat');
summary_csv = fullfile(dataDir, 'figure5a_site_summary.csv');

outBase = fullfile(outputDir, 'Figure5A_SDT_array_bubble_map');

% Bathymetry settings
bathyMin  = -1500;
bathyStep = 100;

% =========================
% Load GEBCO bathymetry
% =========================
lat = double(ncread(ncfile, 'lat'));
lon = double(ncread(ncfile, 'lon'));
[mlon, mlat] = meshgrid(lon, lat);

bathymetry = double(ncread(ncfile, 'elevation'));

% Clip deep bathymetry and put land in separate bin
bathymetry(bathymetry < bathyMin) = bathyMin;
bathymetry(bathymetry > 0) = 1;

% =========================
% Bathymetry colormap
% =========================
waterLvls = bathyMin:bathyStep:0;
lvls = [waterLvls, 1];
nWaterIntervals = numel(waterLvls) - 1;

% Blue bathymetry colormap with pale shallow bin + white land
waterCol = parula(max(nWaterIntervals, 2));
waterCol = flipud(waterCol);
waterCol(end,:) = [0.90 0.95 0.98];
bathyCol = [waterCol; 1 1 1];

% =========================
% Read site coordinates
% =========================
Tsites_all = readtable(sites_csv, 'TextType', 'string');

isSDT = startsWith(upper(strtrim(Tsites_all.Site)), "SDT_");
isSOCALT = upper(strtrim(Tsites_all.Site)) == "SOUTH_T" | upper(strtrim(Tsites_all.Site)) == "SOCAL_T";

Tsites = Tsites_all(isSDT | isSOCALT, :);

% Drop any WWE site if present
isWWE = startsWith(upper(strtrim(Tsites.Site)), "SDT_WWE");
Tsites = Tsites(~isWWE, :);

Tsites.lat_dd = arrayfun(@parseDegMin, Tsites.Lat);
Tsites.lon_dd = arrayfun(@parseDegMin, Tsites.Long);
Tsites.Site_base = stripSDTSuffix(Tsites.Site);

% =========================
% Load UF440 site summary stats
% =========================
summaryTbl = loadSummary(summary_mat, summary_csv);
summaryTbl.Site_base = stripSDTSuffix(summaryTbl.Site);

Presence_percent = nan(height(Tsites),1);
Median_SPL = nan(height(Tsites),1);

for i = 1:height(Tsites)
    sbase = string(Tsites.Site_base(i));
    idx = string(summaryTbl.Site_base) == sbase;
    if any(idx)
        j = find(idx, 1, 'first');
        Presence_percent(i) = summaryTbl.Presence_percent(j);
        Median_SPL(i) = summaryTbl.Median_UF440_sunset_SPL(j);
    end
end

Presence_round = round(Presence_percent);

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
% SPL binning
% =========================
spl = Median_SPL;

nbins = 6;
bin_edges = 71.5:0.5:(71.5 + nbins * 0.5);
bin_labs = bin_edges(1:end-1);

spl_clip = spl;
spl_clip = max(spl_clip, bin_edges(1));
spl_clip = min(spl_clip, bin_edges(end) - eps);

bin_idx = discretize(spl_clip, bin_edges);

% Discrete SPL colormap
cmap7 = jet(7);
cmap6 = cmap7([1 3:7], :);

% =========================
% Figure
% =========================
W_in = 7.5;
H_in = 3.7;

fig = figure(33); clf
set(fig, 'Color', 'w', 'Units', 'inches', 'Position', [1 1 W_in H_in]);
set(fig, 'PaperUnits', 'inches', 'PaperSize', [W_in H_in], ...
         'PaperPosition', [0 0 W_in H_in], 'PaperPositionMode', 'manual');

axPos = [0.08 0.12 0.79 0.82];

% Bathymetry axis
ax1 = axes(fig, 'Position', axPos);
[~, hc] = contourf(ax1, mlon, mlat, bathymetry.', lvls);
set(ax1, 'YDir', 'normal');
hc.LineColor = [0.2 0.2 0.2];
hc.LineWidth = 0.4;
colormap(ax1, bathyCol);
caxis(ax1, [bathyMin 1]);
hold(ax1, 'on');

xlim(ax1, [min(lon) max(lon)]);
ylim(ax1, [min(lat) max(lat)]);
xlabel(ax1, '');
ylabel(ax1, '');
box(ax1, 'on');

% Depth colorbar
cbDepth = colorbar(ax1, 'eastoutside');
cbDepth.TickDirection = 'out';
cbDepth.Label.String = 'Depth (m)';
cbDepth.Ticks = bathyMin:500:0;
cbDepth.TickLabels = compose('%d', bathyMin:500:0);
cbDepth.Position = [0.84 0.20 0.020 0.60];

% Bubble overlay axis
ax2 = axes(fig, 'Position', axPos);
set(ax2, 'Color', 'none', ...
         'XLim', get(ax1, 'XLim'), ...
         'YLim', get(ax1, 'YLim'), ...
         'XTick', [], 'YTick', [], ...
         'Box', 'off');
hold(ax2, 'on');

scatter(ax2, Tsites.lon_dd, Tsites.lat_dd, bubbleA, bin_idx, ...
    'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.0);

colormap(ax2, cmap6);
caxis(ax2, [0.5 nbins + 0.5]);

% SPL colorbar
cbSPL = colorbar(ax2, 'eastoutside');
cbSPL.TickDirection = 'out';
cbSPL.Label.String = 'Median UF440 sunset SPL (dB re 1 \muPa)';
cbSPL.Ticks = 1:nbins;
cbSPL.TickLabels = compose('%.1f', bin_labs);
cbSPL.Position = [0.90 0.20 0.020 0.60];

% =========================
% Save outputs
% =========================
pngFile = strcat(outBase, '.png');
svgFile = strcat(outBase, '.svg');

exportgraphics(fig, pngFile, 'Resolution', 300);
set(fig, 'Renderer', 'painters');
print(fig, svgFile, '-dsvg', '-painters');

fprintf('\nSaved figure:\n  %s\n  %s\n', pngFile, svgFile);

% =========================
% Helper functions
% =========================
function out = stripSDTSuffix(siteStr)
    s = upper(strtrim(string(siteStr)));
    isSDT = startsWith(s, "SDT_");
    out = s;
    out(isSDT) = regexprep(out(isSDT), "_\d+$", "");
end

function dd = parseDegMin(s)
    s = string(s);
    s = strtrim(s);
    parts = split(s);
    hemi  = upper(parts(end));
    dm    = parts(1);
    dmParts = split(dm, "-");
    deg     = str2double(dmParts(1));
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

    req = ["Site", "Presence_percent", "Median_UF440_sunset_SPL"];
    for r = req
        if ~ismember(r, string(summaryTbl.Properties.VariableNames))
            error('Summary file missing required column: %s', r);
        end
    end
    summaryTbl.Site = string(summaryTbl.Site);
end
