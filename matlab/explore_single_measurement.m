%% Single-measurement acquisition check -- fluorescence + reflectance
% Picks ONE *fluo.mat file (its paired *diffuseReflectances.mat is derived
% automatically, same convention as the rest of the pipeline: 'fluo' ->
% 'diffuseReflectances' in the filename) and plots every individual
% acquisition, raw, in its own figure window:
%   1. Reflectance  -- Short,  every acquisition
%   2. Reflectance  -- Long,   every acquisition
%   3. Fluorescence -- 375 nm, every acquisition
%   4. Fluorescence -- 405 nm, every acquisition
%
% Purpose: acquisitions within a single measurement are supposed to be
% repeats of the same spot -- if they visibly drift apart, that usually
% means the probe moved (or something else changed) between acquisitions.
% No normalization or calibration is applied here on purpose, only the raw
% signal, so nothing hides that drift. Each figure can be copied to the
% clipboard ('c') or exported as PNG+PDF ('s') via enable_figure_export.
close all
clear all
projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'matlab', 'functions'))

%% ── PARAMETERS ────────────────────────────────────────────────────────────
dataDir     = 'D:\Lyon thèse\data\Manip cochons Lyon\2024_07_17\System_008_fluo';
min_wl_fluo = 480;
max_wl_fluo = 645;

%% ── SELECTION & LOADING ──────────────────────────────────────────────────
[data, path] = uigetfile('*fluo.mat', 'Select one *fluo.mat measurement', dataDir);
if isequal(data, 0)
    return
end
fluo_path = fullfile(path, data);
label = char(erase(strrep(string(data), "_", " "), ".mat"));

%% ── REFLECTANCE (Short & Long) ────────────────────────────────────────────
reflectance_file = strrep(data, 'fluo', 'diffuseReflectances');
reflectance_path = fullfile(path, reflectance_file);
if ~isfile(reflectance_path)
    warning('No matching reflectance file found (%s) -- skipping the reflectance figures.', reflectance_file);
else
    [~, ~, lambda_refl, ~, ~, ~, WLshort_all, WLlong_all] = load_reflectance(reflectance_path, min_wl_fluo, max_wl_fluo);

    plot_all_acquisitions(lambda_refl, WLshort_all, min_wl_fluo, max_wl_fluo, 'Raw signal (counts)', ...
        [label ' - Short - all acquisitions'], ['single_' label '_refl_short']);
    plot_all_acquisitions(lambda_refl, WLlong_all, min_wl_fluo, max_wl_fluo, 'Raw signal (counts)', ...
        [label ' - Long - all acquisitions'], ['single_' label '_refl_long']);
end

%% ── FLUORESCENCE (375 & 405 nm) ───────────────────────────────────────────
[~, ~, lambda_fluo, fluo375_all, fluo405_all] = load_fluo(fluo_path, min_wl_fluo, max_wl_fluo);

plot_all_acquisitions(lambda_fluo, fluo375_all, min_wl_fluo, max_wl_fluo, 'F.I. (a.u.)', ...
    [label ' - 375 nm - all acquisitions'], ['single_' label '_fluo_375']);
plot_all_acquisitions(lambda_fluo, fluo405_all, min_wl_fluo, max_wl_fluo, 'F.I. (a.u.)', ...
    [label ' - 405 nm - all acquisitions'], ['single_' label '_fluo_405']);

%% ── Local functions ───────────────────────────────────────────────────────
function plot_all_acquisitions(lambda, data, wl_min, wl_max, ylab, figtitle, figname)
%PLOT_ALL_ACQUISITIONS Plots every row of `data` (acquisitions x wavelength)
%   in its own color on a new figure.
f = figure('Position', [100, 100, 700, 450]); hold on
n = size(data, 1);
cmap = turbo(max(n, 1));
for i = 1:n
    plot(lambda, data(i,:), 'LineWidth', 1.5, 'Color', cmap(i,:), 'DisplayName', sprintf('Acquisition %d', i));
end
xlim([wl_min wl_max]); grid on
xlabel('Wavelength (nm)'); ylabel(ylab)
title(figtitle)
legend show
set(gca, 'FontSize', 12, 'FontWeight', 'bold')
enable_figure_export(f, figname);
end
