%% Multi-point biomarker evolution -- compare fluorophore fit results across measurements
% Selects several *fluo.mat files (data points) and runs the fluorophore
% decomposition fit (corrected_fluo_LS_WL_1, same as
% fit_fluo_alpha_fitted_LS_WL_1_unique_acq_function.m) on each one, then
% plots how the results evolve across the selected points (in selection
% order, labelled by filename):
%   1. NADH & FAD evolution -- the two most important biomarkers
%   2. Redox ratio evolution -- its own figure (a ratio, not a % of signal,
%      so it never shares an axis with the fraction plots)
%   3. All 6 fluorophores -- 100%-stacked bar chart per point, one subplot
%      per laser (405 nm | 375 nm)
% Each figure can be copied to the clipboard ('c') or exported as PNG+PDF
% ('s') via enable_figure_export.
%
% Note: corrected_fluo_LS_WL_1 always fits with R=1 (no optical correction --
% same behavior as fit_fluo_alpha_fitted_LS_WL_1_unique_acq_function.m) and
% does not support disabling individual fluorophores.
close all
clear all
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'matlab', 'functions'))

%% ── PARAMETERS ────────────────────────────────────────────────────────────
dataDir            = 'D:\Lyon thèse\data\Manip cochons Lyon\2024_07_17\System_008_fluo';
min_wl_reflectance = 480;
max_wl_reflectance = 645;
min_wl_fluo        = 480;
max_wl_fluo        = 645;
isSmallSpectralon  = true;
name_spectralon    = 'spectralon_000_diffuseReflectances';

fad_path  = fullfile(projectRoot, 'data', 'fluorophores', 'flavine_fluo.mat');
nadh_path = fullfile(projectRoot, 'data', 'fluorophores', 'NADH_fluo_4.mat');

% Same fixed identity colors as python/scripts/_plotting.py::COMPONENT_COLORS
component_colors.FAD      = [0.545, 0.184, 0.788];
component_colors.NADH     = [0.831, 0.627, 0.090];
component_colors.FMN      = [0.102, 0.600, 0.439];
component_colors.Lipo     = [0.180, 0.525, 0.871];
component_colors.PpIX_636 = [0.753, 0.224, 0.169];
component_colors.PpIX_620 = [0.902, 0.494, 0.133];
component_labels.FAD      = 'FAD';
component_labels.NADH     = 'NADH';
component_labels.FMN      = 'protein-bound FMN';
component_labels.Lipo     = 'Lipopigments';
component_labels.PpIX_636 = 'PpIX 636 nm';
component_labels.PpIX_620 = 'PpIX 620 nm';
fluorophore_names = fieldnames(component_colors);  % fixed stacking order

color_375 = [0.0000, 0.4470, 0.7410];
color_405 = [0.8500, 0.3250, 0.0980];

%% ── REFERENCE SPECTRA & FILE SELECTION ────────────────────────────────────
[FAD385, FAD405, ~]   = load_fluo(fad_path, min_wl_fluo, max_wl_fluo);
[NADH385, NADH405, ~] = load_fluo(nadh_path, min_wl_fluo, max_wl_fluo);
fluorophores_385 = [FAD385; NADH385];
fluorophores_405 = [FAD405; NADH405];

[data, path] = uigetfile('*fluo.mat', 'Select the *fluo.mat files (data points) to compare', ...
    'MultiSelect', 'on', dataDir);
if isequal(data, 0)
    return
end
data = convertCharsToStrings(data);
allLabels = erase(strrep(data, "_", " "), ".mat");

%% ── FIT EACH POINT ────────────────────────────────────────────────────────
n = length(data);
labels = strings(1, 0);
frac405 = {}; frac375 = {};
redox405 = []; redox375 = [];

for k = 1:n
    label = allLabels(k);
    fprintf('Fitting %s...\n', label);
    try
        [~, ~, ~, ~, fluorophore, lambda_fit] = corrected_fluo_LS_WL_1( ...
            path, data(k), name_spectralon, isSmallSpectralon, ...
            min_wl_reflectance, max_wl_reflectance, min_wl_fluo, max_wl_fluo, ...
            fluorophores_385, fluorophores_405);
    catch ME
        warning('%s: %s -- skipped.', label, ME.message);
        continue
    end

    labels(end+1) = label; %#ok<SAGROW>
    frac405{end+1} = compute_fractions(fluorophore, lambda_fit, '405'); %#ok<SAGROW>
    frac375{end+1} = compute_fractions(fluorophore, lambda_fit, '385'); %#ok<SAGROW>

    % Redox ratio: FAD / (FAD + NADH), same computation as
    % fit_fluo_alpha_fitted_LS_WL_1_unique_acq_function.m
    FAD405_tot  = sum(fluorophore.flavine_405_exp);
    NADH405_tot = sum(fluorophore.NADH_405_exp);
    FAD375_tot  = sum(fluorophore.flavine_385_exp);
    NADH375_tot = sum(fluorophore.NADH_385_exp);
    redox405(end+1) = FAD405_tot / (FAD405_tot + NADH405_tot); %#ok<SAGROW>
    redox375(end+1) = FAD375_tot / (FAD375_tot + NADH375_tot); %#ok<SAGROW>
end

if isempty(labels)
    error('No file could be fitted successfully.');
end

%% ── FIGURES ───────────────────────────────────────────────────────────────
plot_nadh_fad(labels, frac405, frac375, component_colors);
plot_redox(labels, redox405, redox375, color_405, color_375);

f = figure('Position', [100, 100, figwidth_for(labels, 1400), 600]);
ax_375 = subplot(1,2,1);
ax_405 = subplot(1,2,2);
plot_stacked(ax_375, labels, frac375, '375 nm', fluorophore_names, component_colors, component_labels);
plot_stacked(ax_405, labels, frac405, '405 nm', fluorophore_names, component_colors, component_labels);
sgtitle('Fluorophore composition per point', 'FontWeight', 'bold')
enable_figure_export(f, 'biomarkers_composition');

%% ── Local functions ───────────────────────────────────────────────────────
function fractions = compute_fractions(fluorophore, lambda, laser)
%COMPUTE_FRACTIONS Returns a struct of % fraction of total F.I. per
%   fluorophore for the given laser ('385' or '405'), mirroring
%   python/functions/fluo_model.py::compute_fractions. Handles
%   corrected_fluo_LS_WL_1's inconsistent field-name suffixes directly.
field_map.FAD      = ['flavine_' laser '_exp'];
field_map.NADH     = ['NADH_' laser '_exp'];
field_map.FMN      = ['gaussian_' laser];
field_map.Lipo     = ['lipo_' laser];
field_map.PpIX_636 = ['PpIX_636_' laser];
field_map.PpIX_620 = ['PpIX_620_' laser];

names = fieldnames(field_map);
integrals = struct();
total = 0;
for i = 1:numel(names)
    name = names{i};
    field = field_map.(name);
    if isfield(fluorophore, field) && any(fluorophore.(field) ~= 0)
        value = max(trapz(lambda, fluorophore.(field)), 0);
        integrals.(name) = value;
        total = total + value;
    end
end

fractions = struct();
present = fieldnames(integrals);
for i = 1:numel(present)
    if total > 0
        fractions.(present{i}) = integrals.(present{i}) / total * 100;
    else
        fractions.(present{i}) = 0;
    end
end
end

function values = get_series(fracList, name)
%GET_SERIES Extracts one fluorophore's fraction across all points as a row
%   vector (0 where a point's struct doesn't have that field, e.g. disabled).
values = zeros(1, numel(fracList));
for i = 1:numel(fracList)
    s = fracList{i};
    if isfield(s, name)
        values(i) = s.(name);
    end
end
end

function w = figwidth_for(labels, base_width)
%FIGWIDTH_FOR Widens the figure with the number of points so tick labels
%   have room, mirroring python/scripts/explore_biomarkers_evolution.py.
w = max(base_width, 90 * numel(labels));
end

function plot_nadh_fad(labels, frac405, frac375, component_colors)
x = 1:numel(labels);
f = figure('Position', [100, 100, figwidth_for(labels, 900), 600]); hold on
names = {'NADH', 'FAD'};
for i = 1:numel(names)
    name = names{i};
    color = component_colors.(name);
    plot(x, get_series(frac405, name), '-o',  'Color', color, 'LineWidth', 2, 'DisplayName', [name ' - 405 nm']);
    plot(x, get_series(frac375, name), '--o', 'Color', color, 'LineWidth', 2, 'DisplayName', [name ' - 375 nm']);
end
ax = gca;
set(ax, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45)
ylabel('Fraction of total F.I. (%)')
title('NADH & FAD evolution')
legend show
set(ax, 'FontSize', 12, 'FontWeight', 'bold')
ax.XAxis.FontSize = 8;  % full (untruncated) point labels need a smaller font to fit
enable_figure_export(f, 'biomarkers_nadh_fad');
end

function plot_redox(labels, redox405, redox375, color_405, color_375)
x = 1:numel(labels);
f = figure('Position', [100, 100, figwidth_for(labels, 900), 600]); hold on
plot(x, redox405, '-o', 'Color', color_405, 'LineWidth', 2, 'DisplayName', '405 nm');
plot(x, redox375, '-o', 'Color', color_375, 'LineWidth', 2, 'DisplayName', '375 nm');
ax = gca;
set(ax, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45)
ylabel('Redox ratio  (FAD / (FAD + NADH))')
title('Redox ratio evolution')
legend show
set(ax, 'FontSize', 12, 'FontWeight', 'bold')
ax.XAxis.FontSize = 8;  % full (untruncated) point labels need a smaller font to fit
enable_figure_export(f, 'biomarkers_redox');
end

function plot_stacked(ax, labels, fracList, laserTitle, names, component_colors, component_labels)
%PLOT_STACKED 100%-stacked bar chart of every fluorophore's fraction, one
%   bar per point, on the given axes.
x = 1:numel(labels);
Y = zeros(numel(labels), numel(names));
for i = 1:numel(names)
    Y(:, i) = get_series(fracList, names{i})';
end
active = any(Y ~= 0, 1);
Y = Y(:, active);
activeNames = names(active);

hb = bar(ax, x, Y, 0.6, 'stacked', 'EdgeColor', 'none');
for i = 1:numel(activeNames)
    hb(i).FaceColor = component_colors.(activeNames{i});
    hb(i).DisplayName = component_labels.(activeNames{i});
end
set(ax, 'XTick', x, 'XTickLabel', labels, 'XTickLabelRotation', 45)
ylabel(ax, 'Fraction of total F.I. (%)')
ylim(ax, [0 105])
title(ax, laserTitle)
legend(ax, 'show', 'FontSize', 7)
set(ax, 'FontSize', 12, 'FontWeight', 'bold')
ax.XAxis.FontSize = 8;  % full (untruncated) point labels need a smaller font to fit
end
