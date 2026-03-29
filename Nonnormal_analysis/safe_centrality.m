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
% Authors: Y. Wang, A. N. Montanari, A. E. Motter
% Date: 2026-03-29
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