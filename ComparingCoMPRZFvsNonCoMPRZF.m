% MATLAB Simulation for Coordinated Beamforming (CoMP) vs. Non-CoMP Comparison
% Comparing CoMP RZF with Non-CoMP RZF Beamforming

clear;
clc;
close all;

%% System Parameters
num_users = 50; % Total number of users (I)
num_BS = 4; % Total number of base stations (J) (1 LAP, 3 MBS)
num_antennas_per_BS = 8; % Number of transmit antennas per BS (Nt)
total_antennas = num_BS * num_antennas_per_BS; % Total antennas in the cluster

% Simulation area (square)
area_size_m = 5000; % Simulation area size in meters (5 km x 5 km)

% Random Walk Mobility Parameters
v_min = 1; % Minimum speed (m/s)
v_max = 10; % Maximum speed (m/s)
D_max = 10; % Maximum distance per time slot (meters)

% BS Locations [x, y, z] in meters
center_x = area_size_m / 2;
center_y = area_size_m / 2;
mbs_triangle_side_length = 1800;
mbs_distance_to_center = mbs_triangle_side_length / sqrt(3);

bs_locations = [
    center_x, center_y + mbs_distance_to_center, 30; % MBS 1
    center_x - mbs_triangle_side_length/2, center_y - mbs_distance_to_center/2, 30; % MBS 2
    center_x + mbs_triangle_side_length/2, center_y - mbs_distance_to_center/2, 30; % MBS 3
    center_x, center_y, 2000; % LAP-BS
];

% Channel Model Parameters
fc_MHz = 2545; % Carrier frequency (2.545 GHz)
c = 3e8; % Speed of light
lambda = c / (fc_MHz * 1e6); % Wavelength

% Noise Power
noise_power_dBm = -100; % Noise power in dBm/Hz
noise_power_watt = 10^(noise_power_dBm / 10) * 1e-3; % Noise in Watts

% Path Loss Parameters
mbs_height_m = 30; % MBS height
alpha_NLoS_std_dev_dB = 8; % Shadowing std dev
sigma_LAP = 4; % LAP elevation loss factor
rician_K_dB = 10; % Rician K factor
rician_K = 10^(rician_K_dB / 10); % Linear K factor

% Power Constraints
P_j_dBm = 46; % Per-BS transmit power in dBm
P_j_watt = 10^(P_j_dBm / 10) * 1e-3; % In Watts

% Scheduling & Beamforming Parameters
sus_candidate_set_size = 8; % SUS candidates per BS
SINR_threshold_dB = 17; % RZF SINR threshold
SINR_threshold_linear = 10^(SINR_threshold_dB / 10);
c_init = 1; % RZF initial regularization tuning factor
c_max = 100; % RZF max regularization tuning factor

% Max Simultaneous Users
max_simultaneous_users_comp = total_antennas; % CoMP limit (total spatial streams across cluster)
max_simultaneous_users_per_bs_noncomp = num_antennas_per_BS; % Non-CoMP limit (per BS)

num_time_slots_to_simulate =    1000; % Time slots per user count (for fading averaging)

%% Data Storage
% Store all SINRs for CCDF plots
all_sinrs_comp_rzf_overall = [];
all_sinrs_noncomp_rzf_overall = [];

%% Simulation Execution for 50 Users
fprintf('Starting Simulation for %d Users...\n', num_users);

% Initialize user locations
current_num_users = num_users;
current_user_locations_comp = area_size_m * rand(current_num_users, 2);
current_user_locations_noncomp = area_size_m * rand(current_num_users, 2); % Separate for independence

%% CoMP Simulation
fprintf('  CoMP Simulation (%d time slots)...\n', num_time_slots_to_simulate);

for time_slot = 1:num_time_slots_to_simulate
    % Update User Locations (Random Walk Model)
    theta_i_comp = 2 * pi * rand(current_num_users, 1); % Random direction in [0, 2*pi]
    v_i_comp = v_min + (v_max - v_min) * rand(current_num_users, 1); % Random speed in [v_min, v_max]
    for i = 1:current_num_users
        x_new = current_user_locations_comp(i, 1) + (v_i_comp(i) / v_max) * D_max * cos(theta_i_comp(i));
        y_new = current_user_locations_comp(i, 2) + (v_i_comp(i) / v_max) * D_max * sin(theta_i_comp(i));
        % Boundary reflection
        if x_new < 0
            x_new = -x_new;
            theta_i_comp(i) = pi - theta_i_comp(i);
        elseif x_new > area_size_m
            x_new = 2 * area_size_m - x_new;
            theta_i_comp(i) = pi - theta_i_comp(i);
        end
        if y_new < 0
            y_new = -y_new;
            theta_i_comp(i) = -theta_i_comp(i);
        elseif y_new > area_size_m
            y_new = 2 * area_size_m - y_new;
            theta_i_comp(i) = -theta_i_comp(i);
        end
        current_user_locations_comp(i, :) = [x_new, y_new];
    end
    current_user_locations_comp(:, 1) = max(0, min(area_size_m, current_user_locations_comp(:, 1)));
    current_user_locations_comp(:, 2) = max(0, min(area_size_m, current_user_locations_comp(:, 2)));

    % Channel Generation (including Path Loss and Fading)
    channel_h_ts = cell(current_num_users, num_BS);
    for i = 1:current_num_users
        for j = 1:num_BS
            user_pos_3d = [current_user_locations_comp(i, :), 0];
            bs_pos_3d = bs_locations(j, :);
            d_ij = norm(user_pos_3d - bs_pos_3d);
            if d_ij < 1e-3, d_ij = 1e-3; end

            path_loss_dB_val = 0;
            if j <= 3 % MBS
                alpha_NLoS_dB = alpha_NLoS_std_dev_dB * randn(1);
                path_loss_dB_val = 40 * log10(d_ij / 1000) + (-18 * log10(mbs_height_m) + 21 * log10(fc_MHz) + 80 + alpha_NLoS_dB);
            else % LAP
                d_xy = norm(current_user_locations_comp(i, :) - bs_locations(j, 1:2));
                elevation_angle_rad = atan2(bs_locations(j, 3), d_xy);
                path_loss_free_space_dB = 20 * log10((4 * pi * (fc_MHz * 1e6) * d_ij) / c);
                path_loss_additional_dB = (sigma_LAP - 2) * (1 - elevation_angle_rad / (pi / 2))^2;
                path_loss_dB_val = path_loss_free_space_dB + path_loss_additional_dB;
            end
            channel_gain_linear_ts_ij = 10^(-path_loss_dB_val / 10);

            % Small-scale fading
            if j <= 3 % Rayleigh
                channel_h_ts{i, j} = sqrt(channel_gain_linear_ts_ij) * (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
            else % Rician
                h_LoS = ones(num_antennas_per_BS, 1);
                h_NLoS = (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
                channel_h_ts{i, j} = sqrt(channel_gain_linear_ts_ij) * (sqrt(rician_K / (1 + rician_K)) * h_LoS + sqrt(1 / (1 + rician_K)) * h_NLoS);
            end
        end
    end

    % CoMP Scheduling
    candidate_sets_ts = cell(1, num_BS);
    for j = 1:num_BS
        selected_users_for_bs_j = [];
        available_users = 1:current_num_users;
        current_channel_gains_for_bs_j = zeros(1, current_num_users);
        for u_idx = 1:current_num_users
            current_channel_gains_for_bs_j(u_idx) = norm(channel_h_ts{u_idx, j})^2;
        end
        [~, initial_user_idx] = max(current_channel_gains_for_bs_j(available_users));
        initial_user = available_users(initial_user_idx);
        selected_users_for_bs_j = [selected_users_for_bs_j, initial_user];
        available_users(available_users == initial_user) = [];

        while length(selected_users_for_bs_j) < sus_candidate_set_size && ~isempty(available_users)
            min_max_correlation = inf;
            best_candidate_user = -1;
            for u_candidate = available_users
                max_correlation = 0;
                for u_selected = selected_users_for_bs_j
                    norm_u_candidate = norm(channel_h_ts{u_candidate, j});
                    norm_u_selected = norm(channel_h_ts{u_selected, j});
                    if norm_u_candidate > 0 && norm_u_selected > 0
                        correlation = abs(channel_h_ts{u_candidate, j}' * channel_h_ts{u_selected, j}) / (norm_u_candidate * norm_u_selected);
                        max_correlation = max(max_correlation, correlation);
                    else
                        max_correlation = 1;
                    end
                end
                if max_correlation < min_max_correlation
                    min_max_correlation = max_correlation;
                    best_candidate_user = u_candidate;
                end
            end
            if best_candidate_user ~= -1
                selected_users_for_bs_j = [selected_users_for_bs_j, best_candidate_user];
                available_users(available_users == best_candidate_user) = [];
            else
                break;
            end
        end
        candidate_sets_ts{j} = selected_users_for_bs_j;
    end

    % Conflict Resolution
    user_candidate_counts_ts_comp = zeros(1, current_num_users);
    for i = 1:current_num_users
        for j = 1:num_BS
            if ismember(i, candidate_sets_ts{j})
                user_candidate_counts_ts_comp(i) = user_candidate_counts_ts_comp(i) + 1;
            end
        end
    end
    overlap_users_ts_comp = find(user_candidate_counts_ts_comp > 1);
    non_overlap_users_ts_comp = find(user_candidate_counts_ts_comp == 1);
    resolved_assigned_users_mask_ts_comp = zeros(1, current_num_users);
    resolved_user_bs_assignment_ts_comp = zeros(1, current_num_users);

    for u = overlap_users_ts_comp
        bs_considering_u = find(cellfun(@(c) ismember(u, c), candidate_sets_ts));
        max_gain = -inf;
        best_bs_for_u = -1;
        for j = bs_considering_u
            if norm(channel_h_ts{u, j})^2 > max_gain
                max_gain = norm(channel_h_ts{u, j})^2;
                best_bs_for_u = j;
            end
        end
        if best_bs_for_u ~= -1
            resolved_assigned_users_mask_ts_comp(u) = 1;
            resolved_user_bs_assignment_ts_comp(u) = best_bs_for_u;
        end
    end

    for u = non_overlap_users_ts_comp
        j = find(cellfun(@(c) ismember(u, c), candidate_sets_ts), 1);
        resolved_assigned_users_mask_ts_comp(u) = 1;
        resolved_user_bs_assignment_ts_comp(u) = j;
    end

    % Final selection
    initially_assigned_users_ts_comp = find(resolved_assigned_users_mask_ts_comp == 1);
    if length(initially_assigned_users_ts_comp) > max_simultaneous_users_comp
        user_gain_from_assigned_bs_ts_comp = arrayfun(@(u) norm(channel_h_ts{u, resolved_user_bs_assignment_ts_comp(u)})^2, initially_assigned_users_ts_comp);
        [~, sorted_indices] = sort(user_gain_from_assigned_bs_ts_comp, 'descend');
        S_ts_comp = initially_assigned_users_ts_comp(sorted_indices(1:max_simultaneous_users_comp));
    else
        S_ts_comp = initially_assigned_users_ts_comp;
    end
    num_scheduled_users_ts_comp = length(S_ts_comp);
    s_k_comp = zeros(1, current_num_users);
    s_k_comp(S_ts_comp) = 1; % Binary scheduling variable

    % CoMP RZF Beamforming
    if num_scheduled_users_ts_comp > 0
        H_S_ts_comp = zeros(num_scheduled_users_ts_comp, total_antennas);
        for k = 1:num_scheduled_users_ts_comp
            user_idx = S_ts_comp(k);
            H_S_ts_comp(k, :) = cell2mat(cellfun(@(c) c', channel_h_ts(user_idx, :), 'UniformOutput', false));
        end

        [~, sinrs_comp_rzf_ts] = calculate_rzf_sinr_comp(H_S_ts_comp, s_k_comp(S_ts_comp), P_j_watt, num_antennas_per_BS, noise_power_watt, c_init, c_max, SINR_threshold_linear);
        if ~isempty(sinrs_comp_rzf_ts)
            all_sinrs_comp_rzf_overall = [all_sinrs_comp_rzf_overall, sinrs_comp_rzf_ts];
        end
    end
end

%% Non-CoMP Simulation
fprintf('  Non-CoMP Simulation (%d time slots)...\n', num_time_slots_to_simulate);

for time_slot = 1:num_time_slots_to_simulate
    % Update User Locations (Random Walk Model)
    theta_i_noncomp = 2 * pi * rand(current_num_users, 1);
    v_i_noncomp = v_min + (v_max - v_min) * rand(current_num_users, 1);
    for i = 1:current_num_users
        x_new = current_user_locations_noncomp(i, 1) + (v_i_noncomp(i) / v_max) * D_max * cos(theta_i_noncomp(i));
        y_new = current_user_locations_noncomp(i, 2) + (v_i_noncomp(i) / v_max) * D_max * sin(theta_i_noncomp(i));
        if x_new < 0
            x_new = -x_new;
            theta_i_noncomp(i) = pi - theta_i_noncomp(i);
        elseif x_new > area_size_m
            x_new = 2 * area_size_m - x_new;
            theta_i_noncomp(i) = pi - theta_i_noncomp(i);
        end
        if y_new < 0
            y_new = -y_new;
            theta_i_noncomp(i) = -theta_i_noncomp(i);
        elseif y_new > area_size_m
            y_new = 2 * area_size_m - y_new;
            theta_i_noncomp(i) = -theta_i_noncomp(i);
        end
        current_user_locations_noncomp(i, :) = [x_new, y_new];
    end
    current_user_locations_noncomp(:, 1) = max(0, min(area_size_m, current_user_locations_noncomp(:, 1)));
    current_user_locations_noncomp(:, 2) = max(0, min(area_size_m, current_user_locations_noncomp(:, 2)));

    % Channel Generation
    channel_h_ts = cell(current_num_users, num_BS);
    instantaneous_channel_gain_linear_ts = zeros(current_num_users, num_BS);
    for i = 1:current_num_users
        for j = 1:num_BS
            user_pos_3d = [current_user_locations_noncomp(i, :), 0];
            bs_pos_3d = bs_locations(j, :);
            d_ij = norm(user_pos_3d - bs_pos_3d);
            if d_ij < 1e-3, d_ij = 1e-3; end

            path_loss_dB_val = 0;
            if j <= 3 % MBS
                alpha_NLoS_dB = alpha_NLoS_std_dev_dB * randn(1);
                path_loss_dB_val = 40 * log10(d_ij / 1000) + (-18 * log10(mbs_height_m) + 21 * log10(fc_MHz) + 80 + alpha_NLoS_dB);
            else % LAP
                d_xy = norm(current_user_locations_noncomp(i, :) - bs_locations(j, 1:2));
                elevation_angle_rad = atan2(bs_locations(j, 3), d_xy);
                path_loss_free_space_dB = 20 * log10((4 * pi * (fc_MHz * 1e6) * d_ij) / c);
                path_loss_additional_dB = (sigma_LAP - 2) * (1 - elevation_angle_rad / (pi / 2))^2;
                path_loss_dB_val = path_loss_free_space_dB + path_loss_additional_dB;
            end
            channel_gain_linear_ts_ij = 10^(-path_loss_dB_val / 10);
            instantaneous_channel_gain_linear_ts(i, j) = channel_gain_linear_ts_ij;

            if j <= 3 % Rayleigh
                channel_h_ts{i, j} = sqrt(channel_gain_linear_ts_ij) * (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
            else % Rician
                h_LoS = ones(num_antennas_per_BS, 1);
                h_NLoS = (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
                channel_h_ts{i, j} = sqrt(channel_gain_linear_ts_ij) * (sqrt(rician_K / (1 + rician_K)) * h_LoS + sqrt(1 / (1 + rician_K)) * h_NLoS);
            end
        end
    end

    % Non-CoMP Scheduling (Local SUS per BS)
    scheduled_users_per_bs_ts_noncomp = cell(1, num_BS);
    assigned_bs_for_all_scheduled_users_noncomp_ts = zeros(1, current_num_users);
    s_k_noncomp = zeros(1, current_num_users);

    for j = 1:num_BS
        candidate_users_for_bs_j = sus_scheduler(channel_h_ts, j, current_num_users, sus_candidate_set_size);
        if ~isempty(candidate_users_for_bs_j)
            gains = arrayfun(@(u) instantaneous_channel_gain_linear_ts(u, j), candidate_users_for_bs_j);
            [~, sorted_indices] = sort(gains, 'descend');
            num_to_schedule_j = min(length(candidate_users_for_bs_j), max_simultaneous_users_per_bs_noncomp);
            scheduled_users_j_ts = candidate_users_for_bs_j(sorted_indices(1:num_to_schedule_j));
            final_users_for_bs_j = [];
            for u = scheduled_users_j_ts
                if assigned_bs_for_all_scheduled_users_noncomp_ts(u) == 0
                    final_users_for_bs_j = [final_users_for_bs_j, u];
                    assigned_bs_for_all_scheduled_users_noncomp_ts(u) = j;
                    s_k_noncomp(u) = 1; % Binary scheduling variable
                end
            end
            scheduled_users_per_bs_ts_noncomp{j} = final_users_for_bs_j;
        end
    end
    scheduled_users_ts_noncomp_unique = unique(cell2mat(scheduled_users_per_bs_ts_noncomp));

    % Non-CoMP RZF Beamforming
    if ~isempty(scheduled_users_ts_noncomp_unique)
        W_noncomp_rzf_ts = cell(1, num_BS);
        for j = 1:num_BS
            users_in_Sj = scheduled_users_per_bs_ts_noncomp{j};
            if ~isempty(users_in_Sj)
                H_j_ts = cell2mat(cellfun(@(c) c', channel_h_ts(users_in_Sj, j), 'UniformOutput', false));
                [W_j_rzf, ~] = calculate_rzf_sinr_noncomp(H_j_ts, P_j_watt, noise_power_watt, c_init, c_max, SINR_threshold_linear);
                W_noncomp_rzf_ts{j} = W_j_rzf;
            else
                W_noncomp_rzf_ts{j} = zeros(num_antennas_per_BS, 0);
            end
        end

        % Calculate SINR with inter-cell interference
        current_user_sinrs_ts_noncomp_rzf_slot = [];
        for user_idx = scheduled_users_ts_noncomp_unique
            assigned_bs_idx = assigned_bs_for_all_scheduled_users_noncomp_ts(user_idx);
            if assigned_bs_idx > 0
                users_at_assigned_bs = scheduled_users_per_bs_ts_noncomp{assigned_bs_idx};
                user_local_idx = find(users_at_assigned_bs == user_idx, 1);
                if ~isempty(user_local_idx) && ~isempty(W_noncomp_rzf_ts{assigned_bs_idx})
                    w_k = W_noncomp_rzf_ts{assigned_bs_idx}(:, user_local_idx);
                    h_k_assigned_bs = channel_h_ts{user_idx, assigned_bs_idx};
                    desired_signal_power = s_k_noncomp(user_idx) * abs(h_k_assigned_bs' * w_k)^2;
                    interference_power = 0;
                    for m = 1:num_BS
                        if ~isempty(W_noncomp_rzf_ts{m})
                            h_k_interfering_bs = channel_h_ts{user_idx, m};
                            users_at_m = scheduled_users_per_bs_ts_noncomp{m};
                            if m == assigned_bs_idx
                                for l_idx = 1:length(users_at_m)
                                    if users_at_m(l_idx) ~= user_idx
                                        w_l = W_noncomp_rzf_ts{m}(:, l_idx);
                                        interference_power = interference_power + s_k_noncomp(users_at_m(l_idx)) * abs(h_k_interfering_bs' * w_l)^2;
                                    end
                                end
                            else
                                if ~isempty(users_at_m) && all(users_at_m <= length(s_k_noncomp)) && all(users_at_m > 0)
                                    interference_terms = abs(h_k_interfering_bs' * W_noncomp_rzf_ts{m}).^2;
                                    scheduled_mask = s_k_noncomp(users_at_m);
                                    interference_power = interference_power + sum(scheduled_mask .* interference_terms);
                                end
                            end
                        end
                    end
                    denominator = interference_power + noise_power_watt + eps;
                    current_user_sinrs_ts_noncomp_rzf_slot = [current_user_sinrs_ts_noncomp_rzf_slot, desired_signal_power / denominator];
                end
            end
        end
        valid_sinrs_slot = current_user_sinrs_ts_noncomp_rzf_slot(~isnan(current_user_sinrs_ts_noncomp_rzf_slot) & ~isinf(current_user_sinrs_ts_noncomp_rzf_slot) & current_user_sinrs_ts_noncomp_rzf_slot > 0);
        if ~isempty(valid_sinrs_slot)
            all_sinrs_noncomp_rzf_overall = [all_sinrs_noncomp_rzf_overall, valid_sinrs_slot];
        end
    end
end

fprintf('  Simulations complete for %d users.\n', current_num_users);

fprintf('\n--- All Simulations Complete ---\n');

%% Visualizations

% SINR CCDF Comparison (CoMP RZF vs. Non-CoMP RZF)
fprintf('Generating: SINR CCDF (CoMP RZF vs. Non-CoMP RZF)...\n');
figure;
hold on;
if ~isempty(all_sinrs_comp_rzf_overall)
    [f, x] = ecdf(10*log10(all_sinrs_comp_rzf_overall));
    plot(x, 1-f, 'b-', 'LineWidth', 1, 'DisplayName', 'CoMP (RZF with SUS)');
end
if ~isempty(all_sinrs_noncomp_rzf_overall)
    [f, x] = ecdf(10*log10(all_sinrs_noncomp_rzf_overall));
    plot(x, 1-f, 'r--', 'LineWidth', 1, 'DisplayName', 'Non-CoMP (RZF with SUS)');
end
xlabel('SINR (dB)');
ylabel('CCDF: P(SINR > x)');
title('SINR CCDF (CoMP RZF vs. Non-CoMP RZF)');
legend;
grid on;
hold off;

%% Helper Functions

% SUS Scheduler (used by Non-CoMP)
function selected_users = sus_scheduler(channel_h, bs_idx, num_total_users, candidate_set_size)
    selected_users = [];
    available_users = 1:num_total_users;
    if isempty(available_users), return; end
    gains = arrayfun(@(u) norm(channel_h{u, bs_idx})^2, available_users);
    [~, max_gain_idx] = max(gains);
    initial_user = available_users(max_gain_idx);
    selected_users = [selected_users, initial_user];
    available_users(available_users == initial_user) = [];
    while length(selected_users) < candidate_set_size && ~isempty(available_users)
        min_max_corr = inf;
        best_candidate = -1;
        for u_cand = available_users
            max_corr = 0;
            for u_sel = selected_users
                h_cand = channel_h{u_cand, bs_idx};
                h_sel = channel_h{u_sel, bs_idx};
                corr = abs(h_cand' * h_sel) / (norm(h_cand) * norm(h_sel) + eps);
                if corr > max_corr, max_corr = corr; end
            end
            if max_corr < min_max_corr
                min_max_corr = max_corr;
                best_candidate = u_cand;
            end
        end
        if best_candidate ~= -1
            selected_users = [selected_users, best_candidate];
            available_users(available_users == best_candidate) = [];
        else
            break;
        end
    end
end

% RZF SINR Calculation (CoMP)
function [W_normalized, user_sinrs_linear] = calculate_rzf_sinr_comp(H_S, s_k, P_j_watt, num_antennas_per_BS, noise_power_watt, c_init, c_max, SINR_threshold_linear)
    [num_scheduled_users, total_antennas] = size(H_S);
    num_BS = total_antennas / num_antennas_per_BS;
    W_normalized = []; user_sinrs_linear = [];
    if num_scheduled_users == 0, return; end

    current_c = c_init;
    W_normalized_final = [];
    user_sinrs_linear_final = [];

    while current_c <= c_max
        avg_power_per_user_total = (num_BS * P_j_watt) / max(1, num_scheduled_users);
        delta = current_c * noise_power_watt / avg_power_per_user_total;
        try
            matrix_to_invert = (H_S * H_S' + delta * eye(num_scheduled_users));
            W_unnormalized = H_S' * inv(matrix_to_invert);
        catch
            W_unnormalized = H_S' * pinv(H_S * H_S' + delta * eye(num_scheduled_users));
        end

        if any(isnan(W_unnormalized(:))) || any(isinf(W_unnormalized(:))), return; end

        W_normalized_temp = zeros(total_antennas, num_scheduled_users);
        valid_norm = true;
        for j = 1:num_BS
            ant_idx = ((j-1)*num_antennas_per_BS + 1):(j*num_antennas_per_BS);
            W_j_un = W_unnormalized(ant_idx, :);
            power_j = trace(W_j_un * W_j_un');
            if power_j > P_j_watt
                scale = sqrt(P_j_watt / power_j);
                if isnan(scale) || isinf(scale), valid_norm = false; break; end
                W_normalized_temp(ant_idx, :) = W_j_un * scale;
            else
                W_normalized_temp(ant_idx, :) = W_j_un;
            end
        end
        if ~valid_norm, return; end

        current_sinrs = zeros(1, num_scheduled_users);
        for k = 1:num_scheduled_users
            w_k = W_normalized_temp(:, k);
            h_k = H_S(k, :)';
            desired_power = s_k(k) * abs(h_k' * w_k)^2;
            interference_power = sum(s_k([1:k-1, k+1:end]) .* abs(h_k' * W_normalized_temp(:, [1:k-1, k+1:end])).^2);
            current_sinrs(k) = desired_power / (interference_power + noise_power_watt + eps);
        end

        if any(isnan(current_sinrs)) || any(isinf(current_sinrs)) || any(current_sinrs <= 0)
            % Continue to next 'c'
        else
            if min(current_sinrs) >= SINR_threshold_linear
                W_normalized_final = W_normalized_temp;
                user_sinrs_linear_final = current_sinrs;
                break;
            else
                W_normalized_final = W_normalized_temp;
                user_sinrs_linear_final = current_sinrs;
            end
        end
        current_c = current_c + 0.1;
    end
    W_normalized = W_normalized_final;
    user_sinrs_linear = user_sinrs_linear_final;
end

% RZF SINR Calculation (Non-CoMP)
function [W_normalized, user_sinrs_linear_local] = calculate_rzf_sinr_noncomp(H_j, P_j_watt, noise_power_watt, c_init, c_max, SINR_threshold_linear)
    [num_users_in_Sj, num_antennas_per_BS] = size(H_j);
    W_normalized = zeros(num_antennas_per_BS, 0); user_sinrs_linear_local = [];
    if num_users_in_Sj == 0, return; end

    current_c = c_init;
    W_normalized_final = zeros(num_antennas_per_BS, 0);
    user_sinrs_linear_local_final = [];

    while current_c <= c_max
        avg_power_per_user_j = P_j_watt / max(1, num_users_in_Sj);
        delta_j = current_c * noise_power_watt / avg_power_per_user_j;
        try
            W_unnormalized = H_j' * inv(H_j * H_j' + delta_j * eye(num_users_in_Sj));
        catch
            W_unnormalized = H_j' * pinv(H_j * H_j' + delta_j * eye(num_users_in_Sj));
        end

        if any(isnan(W_unnormalized(:))) || any(isinf(W_unnormalized(:))), return; end

        power_j = trace(W_unnormalized * W_unnormalized');
        if power_j > P_j_watt
            scale = sqrt(P_j_watt / power_j);
            if isnan(scale) || isinf(scale), return; end
            W_normalized_temp = W_unnormalized * scale;
        else
            W_normalized_temp = W_unnormalized;
        end

        current_sinrs_local = zeros(1, num_users_in_Sj);
        for k = 1:num_users_in_Sj
            w_k_j = W_normalized_temp(:, k);
            h_k_j = H_j(k, :)';
            desired_power = abs(h_k_j' * w_k_j)^2;
            interference = sum(abs(h_k_j' * W_normalized_temp(:, [1:k-1, k+1:end])).^2);
            current_sinrs_local(k) = desired_power / (interference + noise_power_watt + eps);
        end

        if any(isnan(current_sinrs_local)) || any(isinf(current_sinrs_local)) || any(current_sinrs_local <= 0)
            % Continue
        else
            if min(current_sinrs_local) >= SINR_threshold_linear
                W_normalized_final = W_normalized_temp;
                user_sinrs_linear_local_final = current_sinrs_local;
                break;
            else
                W_normalized_final = W_normalized_temp;
                user_sinrs_linear_local_final = current_sinrs_local;
            end
        end
        current_c = current_c + 0.1;
    end
    W_normalized = W_normalized_final;
    user_sinrs_linear_local = user_sinrs_linear_local_final;
end