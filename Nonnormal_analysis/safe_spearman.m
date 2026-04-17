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