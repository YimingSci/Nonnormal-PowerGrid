%% add_inertia_to_multiple_generators: Update inertia and primary control for multiple generators
%
% Adds inertia to a set of generators and rescales their primary control
% coefficients proportionally to maintain dynamic consistency.
%
% Inputs:
%   - pant_in  : system struct (pantagruel format)
%   - gen_rows : vector of generator indices (rows in pant.gen)
%   - dM_vec   : vector of inertia increments (same size as gen_rows)
%
% Outputs:
%   - pant_out : updated system struct
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


function pant_out = add_inertia_to_multiple_generators(pant_in, gen_rows, dM_vec)
    pant_out = pant_in;
    for kk = 1:numel(gen_rows)
        dM_add = dM_vec(kk);
        if abs(dM_add) <= 0
            continue;
        end

        idx_gen = gen_rows(kk);
        M_old = pant_out.gen_inertia(idx_gen);
        if M_old <= 0
            error('Generator inertia at row %d is non-positive.', idx_gen);
        end

        M_new = M_old + dM_add;
        k = M_new / M_old;

        pant_out.gen_inertia(idx_gen)   = M_new;
        pant_out.gen_prim_ctrl(idx_gen) = pant_out.gen_prim_ctrl(idx_gen) * k;
    end
end