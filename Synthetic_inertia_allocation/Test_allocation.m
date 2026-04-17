%% Iberian Inertia Allocation and Frequency Response Analysis
%
% This script evaluates how targeted inertia allocation affects transient
% frequency dynamics in a European power grid model (`EUR_2025.mat`).
%
% The workflow includes:
%   - Identifying Iberian renewable generators (PV/WT)
%   - Constructing uniform and nonnormality-based allocation strategies
%   - Simulating disturbances using a Radau5 integrator
%   - Quantifying frequency improvements across Iberian buses
%
% Results are summarized in a compact structure (`Plot_pre`) for analysis.
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


clear; clc;

load('EUR_2025.mat');

%% User settings
M_total_list = 2:0.2:10;
dist_percent_list = 0.1:0.02:1;
p_gain_list = 1 ./ dist_percent_list;

topK         = 5;
dM_test      = 0.01;
dM_step      = 0.10;
gain_floor   = 0;
use_parallel = true;

Sb  = 100;
dt  = 1e-3;
Ndt = 5000;
m_save = 10;
radau_order = 14;
radau_tol   = 1e-6;

interp_max_iter = 500;
interp_alpha    = 1.0;

tol_perf_abs = 0.01;
tol_perf_rel = 0.00;

dst_file = 'New_pre.mat';

%% Parallel pool
if use_parallel
    poolobj = gcp('nocreate');
    if isempty(poolobj)
        parpool;
    end
end

%% Baseline preprocessing
pant0 = pant;

pant0.gen_inertia(619:end,1)   = 0.02;
pant0.gen_prim_ctrl(619:end,1) = 0.005;

pant0.gen_inertia(1:618,1)   = pant0.gen_inertia(1:618,1) * 1;
pant0.gen_prim_ctrl(1:618,1) = pant0.gen_prim_ctrl(1:618,1) * 1;

N_bus = size(pant0.bus,1);

bus_country = normalize_country(pant0.bus_country);
iberia_mask = (bus_country=="ES") | (bus_country=="PT");

is_producing = pant0.gen(:,2) > 0;
id_gen  = pant0.gen(is_producing,1);
id_load = setdiff((1:N_bus)', id_gen);

N_gen  = numel(id_gen);
N_load = numel(id_load);

fault_bus = [3335, 3159, 3324, 3036, 3704, 2537];
fault_mag = [355, 917, 550, 23, 34, 37.5];

%% Candidate buses
[SumNonnormality0, ~, ~, ~] = compute_sum_nonnormality_from_pant(pant0, topK);

Sto_rank = build_sto_rank(pant0, SumNonnormality0);

is_IB_renew = (Sto_rank.GenType=="PV" | Sto_rank.GenType=="WT") & ...
              (Sto_rank.Country=="ES" | Sto_rank.Country=="PT");

Sto_prep = Sto_rank(is_IB_renew,:);
IB_renewable_busnum = Sto_prep.BusNumber;
n_IB_renew = numel(IB_renewable_busnum);

if n_IB_renew == 0
    error('No Iberian renewable buses (PV/WT in ES/PT) found.');
end

candidate_gen_row = zeros(n_IB_renew,1);
for ii = 1:n_IB_renew
    idx_gen = find(pant0.gen(:,1) == IB_renewable_busnum(ii), 1, 'first');
    if isempty(idx_gen)
        error('Bus %d is in candidate set but not found in pant0.gen.', IB_renewable_busnum(ii));
    end
    candidate_gen_row(ii) = idx_gen;
end

fprintf('\n===== Baseline preparation summary =====\n');
fprintf('N_bus                        = %d\n', N_bus);
fprintf('N_gen                        = %d\n', N_gen);
fprintf('N_load                       = %d\n', N_load);
fprintf('Iberian renewable candidates = %d\n', n_IB_renew);
fprintf('M_total sweep                = %.2f : %.2f : %.2f\n', ...
    M_total_list(1), M_total_list(2)-M_total_list(1), M_total_list(end));
fprintf('p_gain sweep                 = %.2f : %.2f : %.2f\n', ...
    p_gain_list(1), p_gain_list(2)-p_gain_list(1), p_gain_list(end));

%% Step 1: allocation
nM = numel(M_total_list);

pant_case1_all = cell(nM,1);
pant_case2_all = cell(nM,1);

nstep_all = zeros(nM,1);

fprintf('\n===== Step 1/2: building allocations for each M_total =====\n');

if use_parallel
    parfor mm = 1:nM
        M_total = M_total_list(mm);

        % Case 1: uniform
        pant_case1 = pant0;
        dM_uniform = M_total / n_IB_renew;
        dM_uniform_vec = dM_uniform * ones(n_IB_renew,1);
        pant_case1 = add_inertia_to_multiple_generators(pant_case1, candidate_gen_row, dM_uniform_vec);

        % Case 2: optimized
        pant_case2 = pant0;
        remaining_M = M_total;
        step_count  = 0;

        while remaining_M > 1e-12
            step_count = step_count + 1;
            this_step  = min(dM_step, remaining_M);

            Sum_cur = compute_sum_nonnormality_from_pant(pant_case2, topK);
            country_cur = normalize_country(pant_case2.bus_country);
            is_IB_cur = (country_cur=="ES") | (country_cur=="PT");
            J_cur_IB = sum(Sum_cur(is_IB_cur));

            gain_raw = zeros(n_IB_renew,1);
            for ii = 1:n_IB_renew
                idx_gen = candidate_gen_row(ii);

                pant_test = pant_case2;
                pant_test = add_inertia_to_one_generator(pant_test, idx_gen, dM_test);

                Sum_test = compute_sum_nonnormality_from_pant(pant_test, topK);
                country_test = normalize_country(pant_test.bus_country);
                is_IB_test = (country_test=="ES") | (country_test=="PT");
                J_test_IB = sum(Sum_test(is_IB_test));

                gain_raw(ii) = J_cur_IB - J_test_IB;
            end

            gain_used = max(gain_raw, gain_floor);

            if sum(gain_used) <= 1e-12
                w_step = ones(n_IB_renew,1) / n_IB_renew;
            else
                w_step = gain_used / sum(gain_used);
            end

            dM_alloc_step = this_step * w_step;
            pant_case2 = add_inertia_to_multiple_generators(pant_case2, candidate_gen_row, dM_alloc_step);

            remaining_M = remaining_M - this_step;
        end

        pant_case1_all{mm} = pant_case1;
        pant_case2_all{mm} = pant_case2;
        nstep_all(mm)      = step_count;
    end
else
    for mm = 1:nM
        M_total = M_total_list(mm);

        % Case 1: uniform
        pant_case1 = pant0;
        dM_uniform = M_total / n_IB_renew;
        dM_uniform_vec = dM_uniform * ones(n_IB_renew,1);
        pant_case1 = add_inertia_to_multiple_generators(pant_case1, candidate_gen_row, dM_uniform_vec);

        % Case 2: optimized
        pant_case2 = pant0;
        remaining_M = M_total;
        step_count  = 0;

        while remaining_M > 1e-12
            step_count = step_count + 1;
            this_step  = min(dM_step, remaining_M);

            Sum_cur = compute_sum_nonnormality_from_pant(pant_case2, topK);
            country_cur = normalize_country(pant_case2.bus_country);
            is_IB_cur = (country_cur=="ES") | (country_cur=="PT");
            J_cur_IB = sum(Sum_cur(is_IB_cur));

            gain_raw = zeros(n_IB_renew,1);
            for ii = 1:n_IB_renew
                idx_gen = candidate_gen_row(ii);

                pant_test = pant_case2;
                pant_test = add_inertia_to_one_generator(pant_test, idx_gen, dM_test);

                Sum_test = compute_sum_nonnormality_from_pant(pant_test, topK);
                country_test = normalize_country(pant_test.bus_country);
                is_IB_test = (country_test=="ES") | (country_test=="PT");
                J_test_IB = sum(Sum_test(is_IB_test));

                gain_raw(ii) = J_cur_IB - J_test_IB;
            end

            gain_used = max(gain_raw, gain_floor);

            if sum(gain_used) <= 1e-12
                w_step = ones(n_IB_renew,1) / n_IB_renew;
            else
                w_step = gain_used / sum(gain_used);
            end

            dM_alloc_step = this_step * w_step;
            pant_case2 = add_inertia_to_multiple_generators(pant_case2, candidate_gen_row, dM_alloc_step);

            remaining_M = remaining_M - this_step;
        end

        pant_case1_all{mm} = pant_case1;
        pant_case2_all{mm} = pant_case2;
        nstep_all(mm)      = step_count;
    end
end

fprintf('Allocation cases finished.\n');

%% Step 2: simulation
nP = numel(p_gain_list);

frac_improved_IB = zeros(nM,nP);

fprintf('\n===== Step 2/2: time-domain simulation for each (M_total, p_gain) =====\n');

combo_list = zeros(nM*nP, 2);
kk = 0;
for mm = 1:nM
    for pp = 1:nP
        kk = kk + 1;
        combo_list(kk,:) = [mm, pp];
    end
end

frac_improved_IB_lin = zeros(size(combo_list,1),1);

if use_parallel
    parfor cc = 1:size(combo_list,1)
        mm = combo_list(cc,1);
        pp = combo_list(cc,2);

        M_total = M_total_list(mm);
        p_gain  = p_gain_list(pp);

        fprintf('Running combo %d / %d | M_total = %.2f | p_gain = %.2f\n', ...
            cc, size(combo_list,1), M_total, p_gain);

        pant_case1 = pant_case1_all{mm};
        pant_case2 = pant_case2_all{mm};

        sim1 = run_time_domain_case( ...
            pant_case1, Sb, dt, Ndt, m_save, ...
            fault_bus, fault_mag, p_gain, ...
            interp_max_iter, interp_alpha, ...
            radau_order, radau_tol);

        sim2 = run_time_domain_case( ...
            pant_case2, Sb, dt, Ndt, m_save, ...
            fault_bus, fault_mag, p_gain, ...
            interp_max_iter, interp_alpha, ...
            radau_order, radau_tol);

        x1 = sim1.peak_bus;
        x2 = sim2.peak_bus;

        tol_here = max(tol_perf_abs, tol_perf_rel * x1);
        improved_mask = (x2 < x1 - tol_here);

        frac_improved_IB_lin(cc) = mean(improved_mask(iberia_mask));
    end
else
    for cc = 1:size(combo_list,1)
        mm = combo_list(cc,1);
        pp = combo_list(cc,2);

        M_total = M_total_list(mm);
        p_gain  = p_gain_list(pp);

        fprintf('Running combo %d / %d | M_total = %.2f | p_gain = %.2f\n', ...
            cc, size(combo_list,1), M_total, p_gain);

        pant_case1 = pant_case1_all{mm};
        pant_case2 = pant_case2_all{mm};

        sim1 = run_time_domain_case( ...
            pant_case1, Sb, dt, Ndt, m_save, ...
            fault_bus, fault_mag, p_gain, ...
            interp_max_iter, interp_alpha, ...
            radau_order, radau_tol);

        sim2 = run_time_domain_case( ...
            pant_case2, Sb, dt, Ndt, m_save, ...
            fault_bus, fault_mag, p_gain, ...
            interp_max_iter, interp_alpha, ...
            radau_order, radau_tol);

        x1 = sim1.peak_bus;
        x2 = sim2.peak_bus;

        tol_here = max(tol_perf_abs, tol_perf_rel * x1);
        improved_mask = (x2 < x1 - tol_here);

        frac_improved_IB_lin(cc) = mean(improved_mask(iberia_mask));
    end
end

%% Recover results
for cc = 1:size(combo_list,1)
    mm = combo_list(cc,1);
    pp = combo_list(cc,2);
    frac_improved_IB(mm,pp) = frac_improved_IB_lin(cc);
end

%% Build output
fprintf('\nData prepared. Building compact plotting struct ...\n');

Plot_pre = struct();

Plot_pre.M_total_list = M_total_list;
Plot_pre.p_gain_list  = p_gain_list;

Plot_pre.nM = numel(M_total_list);
Plot_pre.nP = numel(p_gain_list);

Plot_pre.n_IB_renew   = numel(IB_renewable_busnum);
Plot_pre.tol_perf_abs = tol_perf_abs;
Plot_pre.tol_perf_rel = tol_perf_rel;
Plot_pre.nstep_all    = nstep_all;

Plot_pre.frac_improved_IB = frac_improved_IB;

%% Print summary
fprintf('\n===== Compact plotting summary =====\n');
fprintf('nM = %d, nP = %d\n', Plot_pre.nM, Plot_pre.nP);
fprintf('Number of Iberian renewable buses = %d\n', Plot_pre.n_IB_renew);
fprintf('Mean frac_improved_IB = %.6f\n', mean(frac_improved_IB(:), 'omitnan'));
fprintf('Median frac_improved_IB = %.6f\n', median(frac_improved_IB(:), 'omitnan'));

%% Save
save(dst_file, 'Plot_pre', '-v7');

fprintf('\n===== Finished =====\n');
fprintf('Saved compact file: %s\n', dst_file);

