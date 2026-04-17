%% build_sto_rank: Build a ranked table of bus-level nonnormality and attributes
%
% Sorts buses by the input nonnormality measure and assembles a table
% containing bus number, generator type, inertia, primary control,
% load-frequency coefficient, and country label.
%
% Inputs:
%   - pantagruel      : system struct (pantagruel format)
%   - SumNonnormality : bus-level nonnormality values
%
% Outputs:
%   - Sto_rank : ranked table of bus attributes and nonnormality
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
function Sto_rank = build_sto_rank(pantagruel, SumNonnormality)
    N_bus = size(pantagruel.bus,1);
    [vals_sorted, sort_idx] = sort(SumNonnormality, 'ascend');

    bus_num_all = pantagruel.bus(:,1);
    gen_bus_num = pantagruel.gen(:,1);

    type_sorted        = strings(N_bus,1);
    gen_inertia_sorted = zeros(N_bus,1);
    gen_damping_sorted = zeros(N_bus,1);
    bus_damping_sorted = zeros(N_bus,1);
    bus_num_sorted     = zeros(N_bus,1);
    country_sorted     = strings(N_bus,1);

    gen_type_all = string(pantagruel.gen_type);
    bus_country  = normalize_country(pantagruel.bus_country);

    for ii = 1:N_bus
        idx_bus = sort_idx(ii);

        this_bus_num = bus_num_all(idx_bus);
        bus_num_sorted(ii) = this_bus_num;

        bus_damping_sorted(ii) = pantagruel.load_freq_coef(idx_bus);
        country_sorted(ii)     = bus_country(idx_bus);

        idx_gen = find(gen_bus_num == this_bus_num, 1, 'first');
        if ~isempty(idx_gen)
            type_sorted(ii)        = gen_type_all(idx_gen);
            gen_inertia_sorted(ii) = pantagruel.gen_inertia(idx_gen);
            gen_damping_sorted(ii) = pantagruel.gen_prim_ctrl(idx_gen);
        else
            type_sorted(ii)        = "LOAD";
            gen_inertia_sorted(ii) = 0;
            gen_damping_sorted(ii) = 0;
        end
    end

    Sto_rank = table( ...
        vals_sorted, ...
        bus_num_sorted, ...
        type_sorted, ...
        gen_inertia_sorted, ...
        gen_damping_sorted, ...
        bus_damping_sorted, ...
        country_sorted, ...
        'VariableNames', {'Nonnormality','BusNumber','GenType','GenInertia','GenPrimCtrl','LoadFreqCoef','Country'});
end