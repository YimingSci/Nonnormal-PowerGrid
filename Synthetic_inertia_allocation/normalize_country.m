%% normalize_country: Standardize country labels to uppercase string format
%
% Converts input country identifiers to string format and normalizes them
% by trimming whitespace and converting to uppercase.
%
% Inputs:
%   - raw_country : country labels (cell, string, or char array)
%
% Outputs:
%   - bus_country : normalized string array of country labels
%
% Authors: Y. Wang, A. N. Montanari, A. E. Motter
% Date: 2026-03-29
%
function bus_country = normalize_country(raw_country)
    if iscell(raw_country)
        bus_country = string(raw_country);
    elseif isstring(raw_country)
        bus_country = raw_country;
    elseif ischar(raw_country)
        bus_country = string(cellstr(raw_country));
    else
        error('pant.bus_country data type not recognized.');
    end
    bus_country = upper(strtrim(bus_country));
end