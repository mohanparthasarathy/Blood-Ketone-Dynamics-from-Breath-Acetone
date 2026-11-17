%% NIST acetone spectrum with Fisher-selected channels
% Reads the supplied 67-64-1-IR.jdx file and converts its fractional
% transmittance values to absorbance: A = -log10(T).
%
% The plotted curve is the raw NIST-derived absorbance, not the normalized
% or rescaled spectrum used in the simulated optical measurements.

clear;
clc;

jdxFile = '67-64-1-IR.jdx';

selectedChannels = [1210, 1226, 1450, 1630, 1720, 1736, 1750, 1900];

%% Read the supplied JCAMP-DX file

if ~isfile(jdxFile)
    error(['Cannot find %s. Put the .jdx file in the MATLAB current ' ...
           'folder, or set jdxFile to its full path.'], jdxFile);
end

fileText = fileread(jdxFile);

if isempty(regexp(fileText, ...
        '(?im)^##YUNITS\s*=\s*TRANSMITTANCE\s*$', 'once'))
    error('Expected a JCAMP file with transmittance data.');
end

if isempty(regexp(fileText, ...
        '(?im)^##XYDATA\s*=\s*\(X\+\+\(Y\.\.Y\)\)', 'once'))
    error('This script expects the XYDATA format in the supplied NIST file.');
end

xFactor = readHeaderNumber(fileText, 'XFACTOR');
yFactor = readHeaderNumber(fileText, 'YFACTOR');
deltaX = readHeaderNumber(fileText, 'DELTAX');
expectedPoints = readHeaderNumber(fileText, 'NPOINTS');

lines = regexp(fileText, '\r\n|\n|\r', 'split');

wavenumber = [];
transmittance = [];
readingData = false;

for i = 1:numel(lines)
    line = strtrim(lines{i});

    if ~readingData
        if startsWith(upper(line), '##XYDATA=')
            readingData = true;
        end
        continue;
    end

    if startsWith(line, '##')
        break;
    end

    if isempty(line) || startsWith(line, '$$')
        continue;
    end

    values = sscanf(line, '%f').';

    if numel(values) < 2
        error('Could not read JCAMP data at line %d.', i);
    end

    rowWavenumbers = ...
        (values(1) + (0:numel(values)-2) * deltaX) * xFactor;

    rowTransmittance = values(2:end) * yFactor;

    wavenumber = [wavenumber, rowWavenumbers]; %#ok<AGROW>
    transmittance = [transmittance, rowTransmittance]; %#ok<AGROW>
end

if ~readingData || numel(wavenumber) ~= expectedPoints
    error('The number of data points does not match the JCAMP header.');
end

if any(~isfinite(wavenumber)) || any(~isfinite(transmittance))
    error('The JCAMP data contain non-finite values.');
end

if any(transmittance <= 0) || any(transmittance > 1)
    error(['Expected fractional transmittance values between 0 and 1 ' ...
           'for the supplied file.']);
end

if any(abs(diff(wavenumber) - deltaX * xFactor) > 1e-6)
    error('The wavenumber values are not continuous as expected.');
end

%% Convert transmittance to absorbance and select the plotted range

absorbance = -log10(transmittance);

inRange = wavenumber >= 1000 & wavenumber <= 1900;

% Normalize the plotted spectral range so its maximum absorbance is 1.
absorbance = absorbance / max(absorbance(inRange));

if ~any(inRange) || min(wavenumber) > 1000 || max(wavenumber) < 1900
    error('The NIST spectrum does not cover 1000 to 1900 cm^-1.');
end

wnPlot = wavenumber(inRange);
absPlot = absorbance(inRange);

%% Draw the spectrum and selected channels

fig = figure( ...
    'Color', 'w', ...
    'Position', [100, 100, 1050, 500]);

ax = axes(fig);
ax.Toolbar.Visible = 'off';
hold(ax, 'on');

% Red dotted channel markers.
for i = 1:numel(selectedChannels)
    xline(ax, selectedChannels(i), ':', ...
        'Color', [0.80, 0.15, 0.15], ...
        'LineWidth', 1.2, ...
        'HandleVisibility', 'off');
end

% Draw the spectrum after the markers so the black curve remains visible.
plot(ax, wnPlot, absPlot, ...
    'Color', 'k', ...
    'LineWidth', 1.5);

xlim(ax, [1000, 1900]);

yMax = max(absPlot);
ylim(ax, [0, 1.22 * yMax]);

% Stagger labels vertically to make closely spaced channels easier to read.
labelLevels = [1.08, 1.16, 1.08, 1.08, 1.08, 1.16, 1.08, 1.08];

for i = 1:numel(selectedChannels)
    text(ax, selectedChannels(i), labelLevels(i) * yMax, ...
        sprintf('%d', selectedChannels(i)), ...
        'Color', [0.70, 0.10, 0.10], ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 9, ...
        'Rotation', 90);
end

xlabel(ax, 'Wavenumber (cm^{-1})');
ylabel(ax, 'Normalized absorbance');

set(ax, ...
    'FontName', 'Times New Roman', ...
    'FontSize', 12, ...
    'LineWidth', 0.8, ...
    'Box', 'on');

grid(ax, 'on');
ax.GridAlpha = 0.15;
ax.Layer = 'top';

%% Save files

exportgraphics(fig, 'acetone_spectrum_channels.pdf', ...
    'ContentType', 'vector');

exportgraphics(fig, 'acetone_spectrum_channels.png', ...
    'Resolution', 300);

fprintf('Saved acetone_spectrum_channels.pdf and acetone_spectrum_channels.png\n');

%% Local helper function

function value = readHeaderNumber(fileText, fieldName)
    pattern = ['(?im)^##' fieldName '\s*=\s*([^\r\n]+)'];
    token = regexp(fileText, pattern, 'tokens', 'once');

    if isempty(token)
        error('Missing JCAMP header: %s', fieldName);
    end

    value = str2double(strtrim(token{1}));

    if ~isfinite(value)
        error('Invalid value in JCAMP header: %s', fieldName);
    end
end