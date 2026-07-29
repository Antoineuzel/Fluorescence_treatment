%% Exploration of raw fluorescence spectra (375 & 405 nm)
% Selects one or more *fluo.mat files and plots, for each laser, every
% file's mean +/- 1 standard deviation across its acquisitions (raw and
% normalized), plus a combined view with both lasers on the same axes.
% Each figure can be copied to the clipboard ('c') or exported as PNG+PDF
% ('s') via enable_figure_export.
%
% To inspect a single measurement acquisition-by-acquisition instead (e.g.
% to spot probe movement between repeats), use explore_single_measurement.m.
clear all
close all
addpath(fullfile(fileparts(mfilename('fullpath')), 'functions'))

%% ── PARAMETERS ────────────────────────────────────────────────────────────
dataDir      = 'D:\Lyon thèse\data\Manip cochons Lyon\2024_07_17\System_008_fluo';
min_wl_fluo  = 480;
max_wl_fluo  = 645;

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
[data, path] = uigetfile('*fluo.mat', 'Select one or more files', ...
    'MultiSelect', 'on', dataDir);
data = convertCharsToStrings(data);
legende = erase(strrep(data, "_", " "), ".mat");

n = length(data);
fluo375All = cell(1, n);
fluo405All = cell(1, n);
for i = 1:n
    [~, ~, lambda, fluo375All{i}, fluo405All{i}] = load_fluo(fullfile(path, data(i)), min_wl_fluo, max_wl_fluo);
end

%% ── RAW & NORMALIZED (one figure per laser, mean +/- 1 sigma per file) ───
lasers = {'405', fluo405All, '405 nm'; '375', fluo375All, '375 nm'};

for do_normalize = [false true]
    for r = 1:size(lasers, 1)
        key = lasers{r, 1}; dataAll = lasers{r, 2}; laserTitle = lasers{r, 3};
        f = figure('Position', [100, 100, 700, 450]); hold on
        for i = 1:n
            spec = dataAll{i};
            if do_normalize
                spec = spec ./ sum(spec, 2);
            end
            color_idx = mod(i-1, size(matlab_colors,1)) + 1;
            plot_mean_sigma(gca, lambda, spec, matlab_colors(color_idx,:), legende(i));
        end
        xlim([min_wl_fluo max_wl_fluo]); grid on
        xlabel('Wavelength (nm)');
        if do_normalize
            ylabel('Normalized F.I.'); prefix = 'fluo_norm_';
        else
            ylabel('F.I. (a.u.)'); prefix = 'fluo_raw_';
        end
        title([laserTitle ' excitation'])
        legend show
        set(gca, 'FontSize', 14, 'FontWeight', 'bold')
        enable_figure_export(f, [prefix key]);
    end
end

%% ── BOTH LASERS ON THE SAME AXES (per file mean; color = file, style = laser) ─
f = figure('Position', [100, 100, 1300, 500]);
ax_raw  = subplot(1,2,1); hold on
ax_norm = subplot(1,2,2); hold on
for i = 1:n
    color_idx = mod(i-1, size(matlab_colors,1)) + 1;
    color = matlab_colors(color_idx,:);
    mean405 = mean(fluo405All{i}, 1);
    mean375 = mean(fluo375All{i}, 1);
    plot(ax_raw,  lambda, mean405, '-',  'LineWidth', 2.5, 'Color', color, 'DisplayName', legende(i) + " - 405 nm");
    plot(ax_raw,  lambda, mean375, '--', 'LineWidth', 2.5, 'Color', color, 'DisplayName', legende(i) + " - 375 nm");
    plot(ax_norm, lambda, mean405 / sum(mean405), '-',  'LineWidth', 2.5, 'Color', color, 'DisplayName', legende(i) + " - 405 nm");
    plot(ax_norm, lambda, mean375 / sum(mean375), '--', 'LineWidth', 2.5, 'Color', color, 'DisplayName', legende(i) + " - 375 nm");
end
combo_axes = {ax_raw, 'F.I. (a.u.)', 'Raw'; ax_norm, 'Normalized F.I.', 'Normalized'};
for r = 1:size(combo_axes, 1)
    ax = combo_axes{r,1}; ylab = combo_axes{r,2}; subtitle = combo_axes{r,3};
    xlim(ax, [min_wl_fluo max_wl_fluo]); grid(ax, 'on')
    xlabel(ax, 'Wavelength (nm)'); ylabel(ax, ylab)
    title(ax, subtitle)
    legend(ax, 'show', 'FontSize', 7)
    set(ax, 'FontSize', 13, 'FontWeight', 'bold')
end
sgtitle('375 & 405 nm excitations', 'FontWeight', 'bold')
enable_figure_export(f, 'fluo_combined_lasers');

%% ── Local functions ───────────────────────────────────────────────────────
function plot_mean_sigma(ax, lambda, data, color, lbl)
%PLOT_MEAN_SIGMA Plots the mean +/- 1 standard deviation of `data`
%   (acquisitions x wavelength) on ax, in a single color, labelled `lbl`.
%   The sigma band is excluded from the legend to avoid one "+/- 1 sigma"
%   entry per file.
lambda = lambda(:)';
mu = mean(data, 1); mu = mu(:)';
plot(ax, lambda, mu, 'LineWidth', 2.5, 'Color', color, 'DisplayName', lbl);
if size(data, 1) > 1
    sigma = std(data, 0, 1); sigma = sigma(:)';
    fill(ax, [lambda, fliplr(lambda)], [mu + sigma, fliplr(mu - sigma)], color, ...
        'FaceAlpha', 0.2, 'EdgeColor', 'none', 'HandleVisibility', 'off');
end
end
