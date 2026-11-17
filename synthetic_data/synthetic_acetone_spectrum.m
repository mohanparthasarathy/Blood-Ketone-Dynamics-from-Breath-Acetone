function [wn, absorbance] = synthetic_acetone_spectrum(C_exh_ppm, p)
% NIST-derived relative spectral shape; simulated concentration and scale.
% This reader targets the supplied Coblentz NIST file (plain XYDATA).
% Keep 67-64-1-IR.jdx alongside this function. The 8e-5 scale is an
% illustrative simulation assumption, not an absolute NIST calibration.
persistent shape
wn = linspace(1000, 1900, 451);
if isempty(shape)
    file = fullfile(fileparts(mfilename('fullpath')), '67-64-1-IR.jdx');
    if ~isfile(file)
        error('NIST file not found: %s', file);
    end
    text = fileread(file);
    if isempty(regexp(text, '(?im)^##YUNITS\s*=\s*TRANSMITTANCE\s*$', 'once')) || ...
       isempty(regexp(text, '(?im)^##XYDATA\s*=\s*\(X\+\+\(Y\.\.Y\)\)', 'once'))
        error('Unexpected JCAMP format; this reader requires the supplied NIST file.');
    end
    xfac = header_number(text, 'XFACTOR');
    yfac = header_number(text, 'YFACTOR');
    delta = header_number(text, 'DELTAX');
    expected = header_number(text, 'NPOINTS');
    firstX = header_number(text, 'FIRSTX');
    lastX = header_number(text, 'LASTX');
    lines = regexp(text, '\r\n|\n|\r', 'split');
    reading = false; x = []; y = [];
    for i = 1:numel(lines)
        line = strtrim(lines{i});
        if ~reading
            if startsWith(upper(line), '##XYDATA='); reading = true; end
            continue
        end
        if startsWith(line, '##'); break; end
        if isempty(line) || startsWith(line, '$$'); continue; end
        vals = sscanf(line, '%f').';
        if numel(vals) < 2
            error('Unsupported or malformed XYDATA line %d.', i);
        end
        rowX = (vals(1) + (0:numel(vals)-2)*delta)*xfac;
        x = [x, rowX]; %#ok<AGROW>
        y = [y, vals(2:end)*yfac]; %#ok<AGROW>
    end
    if ~reading || numel(x) ~= expected || any(~isfinite(x)) || any(~isfinite(y))
        error('JCAMP point count or data are invalid.');
    end
    if abs(x(1)-firstX*xfac)>1e-6 || abs(x(end)-lastX*xfac)>1e-6 || ...
       any(abs(diff(x)-delta*xfac)>1e-6)
        error('JCAMP x values do not match the header or are not continuous.');
    end
    % The supplied file has fractional transmission values despite YLABEL
    % reading %% Transmission. Check the actual scale before conversion.
    if any(y <= 0) || any(y > 1)
        error('Expected fractional transmittance in (0,1] for the supplied file.');
    end
    if min(x)>1000 || max(x)<1900
        error('NIST spectrum does not cover the simulation range.');
    end
    A = -log10(y);
    range = x >= 1000 & x <= 1900;
    Amax = max(A(range));
    if ~(Amax > 0); error('Spectrum has no positive absorbance in range.'); end
    shape = interp1(x, A/Amax, wn, 'linear');
end
absorbancePerPpm = 8e-5 * (p.scales.optical_path_m/0.10);
absorbance = absorbancePerPpm * C_exh_ppm(:) * shape;
end

function n = header_number(text, field)
token = regexp(text, ['(?im)^##' field '\s*=\s*([^\r\n]+)'], 'tokens', 'once');
if isempty(token); error('Missing JCAMP header %s.', field); end
n = str2double(strtrim(token{1}));
if ~isfinite(n); error('Invalid JCAMP header %s.', field); end
end
