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
% Authors: Y. Wang, A. N. Montanari, A. E. Motter
% Date: 2026-03-29
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