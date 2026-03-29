%% safe_spearman: Compute Spearman correlation with validity checks
%
% Computes Spearman correlation between two vectors after removing
% invalid entries and ensuring sufficient variability.
%
% Inputs:
%   - x         : input vector
%   - y         : input vector
%   - user_mask : optional logical mask
%
% Outputs:
%   - rho : Spearman correlation (NaN if invalid)
%
% Authors: Y. Wang, A. N. Montanari, A. E. Motter
% Date: 2026-03-29
%
function rho = safe_spearman(x, y, user_mask)
    x = x(:);
    y = y(:);

    if nargin < 3 || isempty(user_mask)
        user_mask = true(size(x));
    end

    mask = user_mask(:) & isfinite(x) & isfinite(y);

    if nnz(mask) < 3
        rho = NaN;
        return;
    end

    x_use = x(mask);
    y_use = y(mask);

    if numel(unique(x_use)) < 2 || numel(unique(y_use)) < 2
        rho = NaN;
        return;
    end

    rho = corr(x_use, y_use, 'type', 'Spearman', 'rows', 'complete');
end