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
% Authors: Y. Wang, A. N. Montanari, A. E. Motter
% Date: 2026-03-29
%
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