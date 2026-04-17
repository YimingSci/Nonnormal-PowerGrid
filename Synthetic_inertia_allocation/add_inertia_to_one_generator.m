%% add_inertia_to_one_generator: Update inertia and primary control of one generator
%
% Adds inertia to a specified generator and rescales its primary control
% coefficient proportionally to preserve dynamic consistency.
%
% Inputs:
%   - pant_in : system struct (pantagruel format)
%   - idx_gen : generator index (row in pant.gen)
%   - dM_add  : inertia increment
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

%
function pant_out = add_inertia_to_one_generator(pant_in, idx_gen, dM_add)
    pant_out = pant_in;
    if abs(dM_add) <= 0
        return;
    end

    M_old = pant_out.gen_inertia(idx_gen);
    if M_old <= 0
        error('Generator inertia at row %d is non-positive.', idx_gen);
    end

    M_new = M_old + dM_add;
    k = M_new / M_old;

    pant_out.gen_inertia(idx_gen)   = M_new;
    pant_out.gen_prim_ctrl(idx_gen) = pant_out.gen_prim_ctrl(idx_gen) * k;
end