%% compute_sum_nonnormality_from_pant: Compute bus-level summed nonnormality
%
% Builds the linearized system matrix, computes its symmetric part, and
% evaluates a bus-level nonnormality measure from the leading eigenmodes.
%
% Inputs:
%   - pantagruel : system struct (pantagruel format)
%   - topK       : number of leading modes included in the sum
%
% Outputs:
%   - SumNonnormality : bus-level summed nonnormality values
%   - A_ext           : linearized system matrix
%   - H               : symmetric part of A_ext
%   - lambda_delta    : difference between the dominant eigenvalues of H and A_ext
%
% Authors: Y. Wang, A. N. Montanari, A. E. Motter
% Date: 2026-03-29
%
function [SumNonnormality, A_ext, H, lambda_delta] = compute_sum_nonnormality_from_pant(pantagruel, topK)
    [A_ext, N_bus] = Build_model(pantagruel);
    H = (A_ext + A_ext') / 2;

    lambda_max_H = eigs(H, 1, 'la');
    lambda_max_A = eigs(A_ext, 1, 'largestreal');
    lambda_delta = lambda_max_H - lambda_max_A;

    topK_eff = min(topK, size(H,1));
    [V_all, D_all] = eigs(H, topK_eff, 'la');
    lambda_all = diag(D_all);

    SumNonnormality = zeros(N_bus,1);
    for kk = 1:topK_eff
        v_tmp = V_all(:,kk) / norm(V_all(:,kk));
        v_nn_total = abs(lambda_all(kk)) * lambda_delta * (v_tmp).^2;
        SumNonnormality = SumNonnormality + v_nn_total(1:N_bus);
    end
    SumNonnormality = log10(SumNonnormality);
end