function enable_figure_export(fig, basename, export_dir)
%ENABLE_FIGURE_EXPORT Binds 'c' (clipboard) and 's' (PNG+PDF export) on fig.
%   enable_figure_export(fig, basename) binds on figure fig two
%   keyboard shortcuts -- MATLAB equivalent of
%   python/scripts/_plotting.py::enable_export :
%     'c' -> copies the figure to the clipboard (copy_figure)
%     's' -> exports as PNG (300 dpi) + vector PDF into export_dir
%   export_dir defaults to: an 'exports' subfolder of the current folder.
if nargin < 3 || isempty(export_dir)
    export_dir = fullfile(pwd, 'exports');
end
if ~exist(export_dir, 'dir')
    mkdir(export_dir);
end
set(fig, 'KeyPressFcn', @(src, evt) on_key(src, evt, basename, export_dir));
end

function on_key(fig, evt, basename, export_dir)
switch evt.Key
    case 'c'
        copy_figure(fig);
        fprintf('Figure copied to clipboard.\n');
    case 's'
        png_path = fullfile(export_dir, [basename '.png']);
        pdf_path = fullfile(export_dir, [basename '.pdf']);
        exportgraphics(fig, png_path, 'Resolution', 300);
        exportgraphics(fig, pdf_path, 'ContentType', 'vector');
        fprintf('Figure exported:\n  %s\n  %s\n', png_path, pdf_path);
end
end
