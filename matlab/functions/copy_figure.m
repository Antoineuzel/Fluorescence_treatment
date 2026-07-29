function copy_figure(fig)
%COPY_FIGURE Copies a figure to the clipboard (vector image).
%   copy_figure(fig) copies fig (defaults to gcf) to the clipboard --
%   pastes directly into Word/PowerPoint, without going through a manual
%   copy/paste of the window. Requires MATLAB R2020a+ (copygraphics).
if nargin < 1 || isempty(fig)
    fig = gcf;
end
try
    copygraphics(fig, 'ContentType', 'vector');
catch
    copygraphics(fig, 'ContentType', 'image');
end
end
