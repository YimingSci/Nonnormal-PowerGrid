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
% -------------------------------------------------------------------------
% Copyright (C) 2026  Y. Wang, A. N. Montanari & A. E. Motter
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software Foundation,
% Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
%
% Last modified by Y. Wang on 2026-03-29

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