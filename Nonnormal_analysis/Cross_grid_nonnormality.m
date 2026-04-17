%% Cross-grid nonnormality–metric correlation analysis
%
% This script evaluates the relationship between nodal nonnormality and
% multiple structural and dynamical metrics across several real-world
% power grids (Poland, Texas, South Korea, Great Britain, Australia, Iceland).
%
% For each grid, the workflow includes:
%   - Constructing the linearized system model
%   - Computing nodal nonnormality from leading symmetric modes
%   - Evaluating Spearman correlations with selected metrics
%   - Aggregating results into a metric-by-grid matrix
%
% The script also:
%   - Generates a polar (radar) plot of |rho| values
%
% Outputs in workspace:
%   - Results        : detailed per-grid results
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
rng(1);

%% =========================================================
%  SETTINGS
%% =========================================================
grid_list = { ...
    struct('name','Poland',        'var_name','case3120sp',                            'ref_freq',50, 'kk_list_max',150), ...
    struct('name','Texas',         'var_name','Texas2k_series24_case3_2024summerpeak', 'ref_freq',60, 'kk_list_max',150), ...
    struct('name','South Korea',   'var_name','KPG193_ver1_5',                         'ref_freq',60, 'kk_list_max',150), ...
    struct('name','Great Britain', 'var_name','GBnetwork',                             'ref_freq',50, 'kk_list_max',150), ...
    struct('name','Australia',     'var_name','snemSA',                                'ref_freq',50, 'kk_list_max',500), ...
    struct('name','Iceland',       'var_name','iceland',                               'ref_freq',50, 'kk_list_max',150) ...
    };

metric_names = { ...
    'Degree', ...
    'Betweenness', ...
    'Closeness', ...
    'Eigenvector', ...
    'PageRank', ...
    'Inertia', ...
    'Damping', ...
    'AbsPinj', ...
    'AvgXDist', ...
    'AbsFiedler'};

metric_labels_short = { ...
    'Deg', 'Bet', 'Close', 'EigC', 'PR', ...
    'H', 'D', '|P|', 'xDist', '|Fied|'};

n_grids   = numel(grid_list);
n_metrics = numel(metric_names);

Results = cell(n_grids,1);

%% =========================================================
%  RUN ALL GRIDS
%% =========================================================
for gg = 1:n_grids
    fprintf('\n====================================================\n');
    fprintf('Running grid %d / %d : %s\n', gg, n_grids, grid_list{gg}.name);
    fprintf('====================================================\n');

    try
        Results{gg} = run_one_grid(grid_list{gg}, metric_names);
    catch ME
        fprintf(2, 'Grid "%s" failed:\n%s\n', grid_list{gg}.name, ME.message);

        tmp = struct();
        tmp.name = grid_list{gg}.name;
        tmp.ok = false;
        tmp.errormsg = ME.message;
        tmp.rho_vec = nan(n_metrics,1);
        Results{gg} = tmp;
    end
end

%% =========================================================
%  OUTPUT MATRIX
%  Row = metric, Column = grid
%% =========================================================
Sto_out = nan(n_metrics, n_grids);
grid_names = cell(1, n_grids);

for gg = 1:n_grids
    grid_names{gg} = Results{gg}.name;
    if isfield(Results{gg}, 'ok') && Results{gg}.ok
        Sto_out(:, gg) = Results{gg}.rho_vec(:);
    end
end

assignin('base', 'Sto_out', Sto_out);
assignin('base', 'metric_names', metric_names);
assignin('base', 'metric_labels_short', metric_labels_short);
assignin('base', 'grid_names', grid_names);
assignin('base', 'Results', Results);

fprintf('\n\n==================== GLOBAL SUMMARY ====================\n');
T_summary = array2table(Sto_out, ...
    'VariableNames', matlab.lang.makeValidName(grid_names), ...
    'RowNames', metric_names);
disp(T_summary);

%% =========================================================
%  FIGURE 20: RADAR / POLAR MAP
%  Radius = |rho|
%% =========================================================
figure(20); clf;
set(gcf, 'Color', 'w', 'Position', [150, 80, 860, 760]);

theta_deg_full = linspace(0, 360, n_metrics + 1);
theta_full = deg2rad(theta_deg_full);

pax = polaraxes;
hold(pax, 'on');
pax.ThetaZeroLocation = 'top';
pax.ThetaDir = 'clockwise';
pax.RLim = [0 1];
pax.RTick = [0.25 0.5 0.75 1.0];
pax.RTickLabel = {'0.25','0.5','0.75','1.0'};
pax.ThetaTick = theta_deg_full(1:end-1);
pax.ThetaTickLabel = metric_labels_short;
pax.FontSize = 10;
pax.LineWidth = 1.0;

title('Figure 20. Polar map of |Spearman \rho| across selected metrics', ...
    'FontWeight','bold', 'FontSize', 13);

for gg = 1:n_grids
    rho_abs = abs(Sto_out(:, gg));
    rho_plot = [rho_abs(:); rho_abs(1)];

    polarplot(theta_full, rho_plot, '-o', ...
        'LineWidth', 1.8, ...
        'MarkerSize', 4.8, ...
        'DisplayName', grid_names{gg});
end

legend('Location','southoutside', 'NumColumns', 3, 'Box','off');

%% =========================================================
%  NEW PART:
%  OUTPUT r/theta MATRICES
%
%  Each Sto_* matrix is 10x2:
%      col 1 = r = |rho|
%      col 2 = theta (degree)
%% =========================================================
theta_metric_deg = linspace(0, 360, n_metrics + 1);
theta_metric_deg = theta_metric_deg(1:end-1).';   % 10x1

for gg = 1:n_grids
    if isfield(Results{gg}, 'ok') && Results{gg}.ok
        r = abs(Results{gg}.rho_vec(:));
        theta = theta_metric_deg;

        Sto_rt = [r, theta];

        switch Results{gg}.name
            case 'Poland'
                var_name_rt = 'Sto_Poland';
            case 'Texas'
                var_name_rt = 'Sto_Texas';
            case 'South Korea'
                var_name_rt = 'Sto_SouthKorea';
            case 'Great Britain'
                var_name_rt = 'Sto_GreatBritain';
            case 'Australia'
                var_name_rt = 'Sto_Australia';
            case 'Iceland'
                var_name_rt = 'Sto_Iceland';
            otherwise
                var_name_rt = ['Sto_' matlab.lang.makeValidName(Results{gg}.name)];
        end

        assignin('base', var_name_rt, Sto_rt);

        fprintf('\n%s has been created in workspace: %s\n', ...
            var_name_rt, mat2str(size(Sto_rt)));
    end
end