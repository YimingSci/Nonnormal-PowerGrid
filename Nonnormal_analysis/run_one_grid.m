%% run_one_grid: Compute nonnormality and metric correlations for one grid
%
% Loads a grid model, constructs the linearized dynamics, computes
% nodal nonnormality, and evaluates Spearman correlations between
% nonnormality and selected structural/dynamic metrics.
%
% Inputs:
%   - grid_cfg     : struct with grid settings and model name
%   - metric_names : list of metric names (for reporting)
%
% Outputs:
%   - R : result struct containing rho values and node-level data
%
% Authors: Y. Wang, A. N. Montanari, A. E. Motter
% Date: 2026-03-29
%
function R = run_one_grid(grid_cfg, metric_names)

    R = struct();
    R.name = grid_cfg.name;
    R.ok = false;

    kk_list_max = grid_cfg.kk_list_max;

    %% --------------------- STEP 0: load grid ---------------------
    try
        mpc_raw = eval(grid_cfg.var_name);
    catch
        error('Cannot evaluate grid case "%s". Please check the model name.', grid_cfg.var_name);
    end
    mpc_raw.ref_freq = grid_cfg.ref_freq;

    %% ===================== STEP 1: dynamic parameters =====================
    [p, ~] = EN_model(mpc_raw);

    ngi_raw = size(mpc_raw.gen,1);
    H_taka  = p.H;
    D_taka  = p.D;
    OmegaR  = p.omega_R;

    if numel(H_taka) ~= ngi_raw
        H_ext = mean(H_taka) * ones(ngi_raw,1);
        D_ext = mean(D_taka) * ones(ngi_raw,1);
    else
        H_ext = H_taka;
        D_ext = D_taka;
    end

    mpc = mpc_raw;
    mpc.gen_inertia    = 2 * H_ext / OmegaR;
    mpc.gen_prim_ctrl  = D_ext;
    mpc.load_freq_coef = mean(mpc.gen_prim_ctrl) * ones(size(mpc.bus,1),1) * 0.1;

    pantagruel = mpc;

    %% ===================== STEP 2: build model =====================
    [~, A_ext, N_bus] = Build_model(pantagruel);
    Hmat = (A_ext + A_ext') / 2;

    lambda_max_H = eigs(Hmat, 1, 'la');
    lambda_max_A = eigs(A_ext, 1, 'largestreal');
    lambda_delta = lambda_max_H - lambda_max_A;

    %% ===================== STEP 3: nodal nonnormality =====================
    kk_list = 1:kk_list_max;
    [V_all, D_all] = eigs(Hmat, max(kk_list), 'la');
    lambda_all = diag(D_all);

    SumReactivity_raw = zeros(N_bus,1);

    for kk = kk_list
        v_tmp = V_all(:, kk) / norm(V_all(:, kk));
        v_re_total = abs(lambda_all(kk)) * lambda_delta * v_tmp;

        % first N_bus states = frequency subspace
        freq_bus = v_re_total(1:N_bus);
        SumReactivity_raw = SumReactivity_raw + abs(freq_bus).^2;
    end

    SumReactivity_log = log10(SumReactivity_raw + 1e-300);
    SumReactivity_norm = (SumReactivity_log - min(SumReactivity_log)) / ...
                         (max(SumReactivity_log) - min(SumReactivity_log) + 1e-12);

    %% ===================== STEP 4: bus / graph =====================
    define_constants
    bus_id_ext = pantagruel.bus(:, BUS_I);

    [tf1, fromBus] = ismember(pantagruel.branch(:, F_BUS), bus_id_ext);
    [tf2, toBus]   = ismember(pantagruel.branch(:, T_BUS), bus_id_ext);
    if any(~tf1) || any(~tf2)
        error('Branch bus id not found in BUS_I.');
    end

    x_line = abs(pantagruel.branch(:, BR_X));
    x_line(x_line < 1e-6) = 1e-6;

    G_unw = graph(fromBus, toBus);
    G_x   = graph(fromBus, toBus, x_line);   % weighted by reactance distance

    %% ===================== STEP 5: structural metrics =====================
    degree_bus = accumarray([fromBus; toBus], 1, [N_bus, 1], @sum, 0);
    betweenness_bus = safe_centrality(G_unw, 'betweenness', N_bus);
    closeness_bus   = safe_centrality(G_unw, 'closeness',   N_bus);
    eigencent_bus   = safe_centrality(G_unw, 'eigenvector', N_bus);
    pagerank_bus    = safe_centrality(G_unw, 'pagerank',    N_bus);

    %% ===================== STEP 6: inertia / damping / injections =====================
    gen_bus_ids = pantagruel.gen(:, GEN_BUS);
    gen_inertia_each = pantagruel.gen_inertia(:);
    gen_ctrl_each    = pantagruel.gen_prim_ctrl(:);

    [tf_gen, gen_bus_pos] = ismember(gen_bus_ids, bus_id_ext);
    if any(~tf_gen)
        error('Generator bus id not found in BUS_I.');
    end

    inertia_bus = accumarray(gen_bus_pos, gen_inertia_each, [N_bus, 1], @sum, 0);

    gen_ctrl_bus = accumarray(gen_bus_pos, gen_ctrl_each, [N_bus, 1], @sum, 0);
    gen_bus_mask = inertia_bus > 0;

    damping_bus = gen_ctrl_bus + pantagruel.load_freq_coef(:);

    Pd = pantagruel.bus(:, PD);
    Pg_bus = accumarray(gen_bus_pos, pantagruel.gen(:, PG), [N_bus,1], @sum, 0);
    Pinj_abs = abs(Pg_bus - Pd);

    %% ===================== STEP 7: electrical distance =====================
    try
        D_x = distances(G_x);
        D_x(~isfinite(D_x)) = NaN;
        D_x(1:N_bus+1:end) = NaN;
        AvgXDist_bus = nanmean(D_x, 2);
    catch
        warning('Weighted shortest-path metric failed for %s.', grid_cfg.name);
        AvgXDist_bus = nan(N_bus,1);
    end

    %% ===================== STEP 8: weighted Laplacian spectral metrics =====================
    W = sparse([fromBus; toBus], [toBus; fromBus], [1./x_line; 1./x_line], N_bus, N_bus);
    dW = sum(W, 2);
    Lw = spdiags(dW, 0, N_bus, N_bus) - W;

    try
        k_lap = min([30, N_bus-1]);
        [Ulap, Dlap] = eigs(Lw, k_lap + 1, 'smallestabs');
        lam = real(diag(Dlap));

        [lam, idx_sort] = sort(lam, 'ascend');
        Ulap = Ulap(:, idx_sort);

        tol_lam = 1e-8;
        keep = lam > tol_lam;
        lam_nz = lam(keep);
        U_nz   = Ulap(:, keep);

        if numel(lam_nz) >= 1
            AbsFiedler_bus = abs(U_nz(:,1));
        else
            AbsFiedler_bus = nan(N_bus,1);
        end
    catch
        warning('Weighted Laplacian spectral metric failed for %s.', grid_cfg.name);
        AbsFiedler_bus = nan(N_bus,1);
    end

    %% ===================== STEP 9: compute all rhos =====================
    rho_degree       = safe_spearman(degree_bus,      SumReactivity_norm, true(N_bus,1));
    rho_betweenness  = safe_spearman(betweenness_bus, SumReactivity_norm, true(N_bus,1));
    rho_closeness    = safe_spearman(closeness_bus,   SumReactivity_norm, true(N_bus,1));
    rho_eigencent    = safe_spearman(eigencent_bus,   SumReactivity_norm, true(N_bus,1));
    rho_pagerank     = safe_spearman(pagerank_bus,    SumReactivity_norm, true(N_bus,1));

    rho_inertia      = safe_spearman(inertia_bus,     SumReactivity_norm, gen_bus_mask);
    rho_damping      = safe_spearman(damping_bus,     SumReactivity_norm, true(N_bus,1));
    rho_Pinj_abs     = safe_spearman(Pinj_abs,        SumReactivity_norm, true(N_bus,1));
    rho_AvgXDist     = safe_spearman(AvgXDist_bus,    SumReactivity_norm, true(N_bus,1));
    rho_AbsFiedler   = safe_spearman(AbsFiedler_bus,  SumReactivity_norm, true(N_bus,1));

    rho_vec = [ ...
        rho_degree; ...
        rho_betweenness; ...
        rho_closeness; ...
        rho_eigencent; ...
        rho_pagerank; ...
        rho_inertia; ...
        rho_damping; ...
        rho_Pinj_abs; ...
        rho_AvgXDist; ...
        rho_AbsFiedler];

    %% ===================== STEP 10: summary =====================
    fprintf('\n========== SUMMARY: %s ==========\n', grid_cfg.name);
    fprintf('N_bus = %d\n', N_bus);
    for ii = 1:numel(metric_names)
        fprintf('%-15s : %+8.4f\n', metric_names{ii}, rho_vec(ii));
    end
    fprintf('=============================================\n');

    %% ===================== save =====================
    R.ok = true;
    R.N_bus = N_bus;
    R.SumReactivity_log  = SumReactivity_log;
    R.SumReactivity_norm = SumReactivity_norm;
    R.rho_vec = rho_vec;

    R.degree_bus       = degree_bus;
    R.betweenness_bus  = betweenness_bus;
    R.closeness_bus    = closeness_bus;
    R.eigencent_bus    = eigencent_bus;
    R.pagerank_bus     = pagerank_bus;
    R.inertia_bus      = inertia_bus;
    R.damping_bus      = damping_bus;
    R.Pinj_abs         = Pinj_abs;
    R.AvgXDist_bus     = AvgXDist_bus;
    R.AbsFiedler_bus   = AbsFiedler_bus;

    R.G_unw = G_unw;
    R.G_x   = G_x;
    R.fromBus = fromBus;
    R.toBus   = toBus;
end