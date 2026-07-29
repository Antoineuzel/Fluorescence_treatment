%% Exploration of diffuse reflectance (short/long wavelength)
% Selects one or more *diffuseReflectances.mat files and plots four stages,
% in order:
%   1. Raw signal        -- straight from the instrument, no spectralon division yet
%   2. Reflectance        -- raw signal / spectralon source spectrum
%   3. Normalized          -- stage 2, each spectrum normalized to area = 1
%   4. Scaled reflectance  -- raw signal / (spectralon * scale factor), in cm^-2
% Stages 2-4 require the spectralon + scale calibration files; if missing,
% only the raw signal is shown. Each figure can be copied to the clipboard
% ('c') or exported as PNG+PDF ('s') via enable_figure_export.
close all
clear all
warning('off', 'all')
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'matlab', 'functions'))

%% ── PARAMETERS ────────────────────────────────────────────────────────────
dataDir             = 'D:\Lyon thèse\data\Manip cochons Lyon\mixed_pig';
min_wl              = 480;
max_wl              = 645;
isSmallSpectralon   = true;
calibrate           = true;   % false -> only the raw signal is shown
name_spectralon     = 'spectralon_000_diffuseReflectances';  % default name if the .mat file does not specify one

file_spectralon_theo = fullfile(projectRoot, 'data', 'spectralon', 'Reflectance_values_array.txt');
scale_path            = fullfile(projectRoot, 'data', 'probe_calibration', 'scale_new_theory_intralipids.mat');

matlab_colors = [
    0.0000, 0.4470, 0.7410;  % blue
    0.8500, 0.3250, 0.0980;  % orange
    0.9290, 0.6940, 0.1250;  % yellow
    0.4940, 0.1840, 0.5560;  % purple
    0.4660, 0.6740, 0.1880;  % green
    0.3010, 0.7450, 0.9330;  % cyan
    0.6350, 0.0780, 0.1840;  % dark red
    0.8500, 0.5000, 0.5000   % pink
];

%% ── SELECTION & LOADING ──────────────────────────────────────────────────
[files, path] = uigetfile('*diffuseReflectances.mat', 'Select one or more files', ...
    'MultiSelect', 'on', dataDir);
files = convertCharsToStrings(files);
legende = erase(strrep(files, "_", " "), "diffuseReflectances.mat");

can_calibrate = calibrate && isfile(file_spectralon_theo) && isfile(scale_path);
if calibrate && ~can_calibrate
    warning('Calibration files not found (%s) -- displaying raw signal only.', ...
        file_spectralon_theo);
end
if can_calibrate
    scale = load(scale_path);
end

n = length(files);
RawShortAll = cell(1, n); RawLongAll = cell(1, n);
ReflShortAll = cell(1, n); ReflLongAll = cell(1, n);
ScaledShortAll = cell(1, n); ScaledLongAll = cell(1, n);
lambda_wl = [];

for k = 1:n
    file_WL = fullfile(path, files(k));
    [WLShort, WLLong, lambda_wl, WLight, idx_min, idx_max] = load_reflectance(file_WL, min_wl, max_wl);
    RawShortAll{k} = WLShort;
    RawLongAll{k}  = WLLong;

    if ~can_calibrate
        continue
    end
    if ~isfield(WLight, 'lambda')
        WLight.lambda = WLight.raw_lambda;
    end

    if isfield(WLight, 'name_spectralon')
        file_spectralon_exp = fullfile(path, WLight.name_spectralon);
    else
        file_spectralon_exp = fullfile(path, [name_spectralon '.mat']);
    end
    if ~isfile(file_spectralon_exp)
        warning('Spectralon not found for %s -- reflectance skipped for this file.', files(k));
        continue
    end

    [SsourceShort, SsourceLong] = calibration_spectralon(file_spectralon_exp, file_spectralon_theo);
    if isSmallSpectralon
        SsourceShort = SsourceShort / 0.8;
        SsourceLong  = SsourceLong / 0.8;
    end
    unscaled_short = SsourceShort(idx_min:idx_max);
    unscaled_long  = SsourceLong(idx_min:idx_max);

    scale_short = interp1(scale.Scale.Short.WL, scale.Scale.Short.Sy, lambda_wl);
    scale_long  = interp1(scale.Scale.Long.WL,  scale.Scale.Long.Sy,  lambda_wl);

    ReflShortAll{k}   = WLShort ./ unscaled_short;
    ReflLongAll{k}    = WLLong  ./ unscaled_long;
    ScaledShortAll{k} = WLShort ./ (unscaled_short .* scale_short);
    ScaledLongAll{k}  = WLLong  ./ (unscaled_long  .* scale_long);
end

%% ── 1. RAW SIGNAL ─────────────────────────────────────────────────────────
plot_short_long(RawShortAll, RawLongAll, lambda_wl, legende, matlab_colors, [min_wl max_wl], ...
    'Raw signal (counts)', 'reflectance_raw_signal');

%% ── 2. REFLECTANCE (unscaled) ─────────────────────────────────────────────
if can_calibrate
    plot_short_long(ReflShortAll, ReflLongAll, lambda_wl, legende, matlab_colors, [min_wl max_wl], ...
        'Reflectance', 'reflectance_unscaled');
end

%% ── 3. NORMALIZED REFLECTANCE ─────────────────────────────────────────────
if can_calibrate
    normShort = cellfun(@(v) v / sum(v), ReflShortAll, 'UniformOutput', false);
    normLong  = cellfun(@(v) v / sum(v), ReflLongAll,  'UniformOutput', false);
    plot_short_long(normShort, normLong, lambda_wl, legende, matlab_colors, [min_wl max_wl], ...
        'Normalized reflectance', 'reflectance_normalized');
end

%% ── 4. SCALED REFLECTANCE (last) ──────────────────────────────────────────
if can_calibrate
    plot_short_long(ScaledShortAll, ScaledLongAll, lambda_wl, legende, matlab_colors, [min_wl max_wl], ...
        'Scaled reflectance (cm^{-2})', 'reflectance_scaled');
end

%% ── Local functions ───────────────────────────────────────────────────────
function plot_short_long(shortAll, longAll, lambda_wl, legende, matlab_colors, xl, ylab, figname)
%PLOT_SHORT_LONG Plots a Short/Long pair of cell arrays (one entry per file)
%   as two subplots, skipping files with empty data (e.g. missing spectralon).
f = figure('Position', [100, 100, 1100, 500]);
ax_short = subplot(1,2,1); hold on
ax_long  = subplot(1,2,2); hold on
for k = 1:numel(shortAll)
    if isempty(shortAll{k})
        continue
    end
    color_idx = mod(k-1, size(matlab_colors,1)) + 1;
    plot(ax_short, lambda_wl, shortAll{k}, 'LineWidth', 3, 'Color', matlab_colors(color_idx,:), 'DisplayName', legende(k));
    plot(ax_long,  lambda_wl, longAll{k},  'LineWidth', 3, 'Color', matlab_colors(color_idx,:), 'DisplayName', legende(k));
end
axes_list = {ax_short, 'Short'; ax_long, 'Long'};
for i = 1:size(axes_list, 1)
    ax = axes_list{i,1};
    xlim(ax, xl); grid(ax, 'on')
    xlabel(ax, 'Wavelength (nm)'); ylabel(ax, ylab)
    title(ax, axes_list{i,2})
    legend(ax, 'show', 'FontSize', 8)
    set(ax, 'FontSize', 13, 'FontWeight', 'bold')
end
enable_figure_export(f, figname);
end
