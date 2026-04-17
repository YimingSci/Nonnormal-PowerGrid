%% safe_centrality: Compute graph centrality with error handling
%
% Computes a specified centrality measure for a graph and returns
% a valid vector of size N_bus. If computation fails, returns NaNs.
%
% Inputs:
%   - G           : graph object
%   - method_name : centrality type (e.g., 'betweenness')
%   - N_bus       : expected number of nodes
%
% Outputs:
%   - out : centrality vector (Nx1)
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
function out = safe_centrality(G, method_name, N_bus)
    try
        out = centrality(G, method_name);
        out = out(:);
        if numel(out) ~= N_bus
            out = nan(N_bus,1);
        end
    catch
        warning('centrality(%s) failed.', method_name);
        out = nan(N_bus,1);
    end
end