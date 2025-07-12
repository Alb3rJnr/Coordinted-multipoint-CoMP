% MATLAB Simulation to compare SINR variation with increasing number of users
% Compares CoMP RZF and Non-CoMP RZF performance


clear;
clc;
close all;

%% System Parameters (Fixed for all user counts)
num_BS = 4; % Total number of base stations (J) (1 LAP, 3 MBS)
num_antennas_per_BS = 8; % Number of transmit antennas per BS (Nt)
total_antennas = num_BS * num_antennas_per_BS; % Total antennas in the cluster

% Simulation area (square)
area_size_m = 5000; % Simulation area size in meters (5 km x 5 km)

% Random Walk Mobility Parameters
v_min = 1; % Minimum speed (m/s)
v_max = 10; % Maximum speed (m/s)
D_max = 10; % Maximum distance per time slot (meters)

% BS Locations [x, y, z] in meters (Fixed)
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

% Channel Model Parameters (Fixed)
fc_MHz = 2545; % Carrier frequency in MHz (2.545 GHz)
c = 3e8; % Speed of light in m/s
lambda = c / (fc_MHz * 1e6); % Wavelength
noise_power_dBm = -100; % Noise power in dBm
noise_power_watt = 10^(noise_power_dBm / 10) * 1e-3; % Noise power in Watts (N0)

% Path Loss Model Parameters (Fixed)
mbs_height_m = 30; % MBS height
alpha_NLoS_std_dev_dB = 6; % Standard deviation for alpha_NLoS (shadowing)
sigma_LAP = 4; % Elevation-dependent loss factor for LAP
rician_K_dB = 10; % Rician K factor in dB for LAP-UE links
rician_K = 10^(rician_K_dB / 10); % Rician K factor (linear)

% Power Constraints (Fixed)
P_j_dBm = 46; % Per-BS transmit power in dBm
P_j_watt = 10^(P_j_dBm / 10) * 1e-3; % Per-BS transmit power in Watts
P_total_watt = num_BS * P_j_watt; % Total transmit power in Watts

% Scheduling & Beamforming Parameters (Fixed)
sus_candidate_set_size = 8; % Number of users each BS initially selects as candidates (CoMP SUS)
SINR_threshold_dB = 17; % Minimum required SINR in dB (gamma_thre) for RZF convergence
SINR_threshold_linear = 10^(SINR_threshold_dB / 10); % Minimum required SINR (linear)
c_init = 1; % Initial regularization parameter scaling factor (RZF)
c_max = 100; % Maximum regularization parameter scaling factor (RZF)

% Maximum number of simultaneously scheduled users across the cluster (for CoMP feasibility)
max_simultaneous_users_comp = total_antennas;

% Maximum number of simultaneously scheduled users *per BS* (for Non-CoMP feasibility)
max_simultaneous_users_per_bs_noncomp = num_antennas_per_BS;

%% Simulation Parameters (Varying Number of Users and Fixed Time Slots)
num_users_range = [10, 20, 30, 40, 50]; % Range of user numbers to simulate
num_time_slots_to_simulate =   1000; % Fixed number of time slots for each user count

%% Data Storage for SINR vs. Number of Users
median_sinr_comp_rzf_vs_users = zeros(1, length(num_users_range));
median_sinr_noncomp_rzf_vs_users = zeros(1, length(num_users_range));

%% Simulation Loop (Iterate over different numbers of users)
fprintf('Starting Simulation across different numbers of users...\n');

for user_count_idx = 1:length(num_users_range)
    current_num_users = num_users_range(user_count_idx);
    fprintf('\n--- Simulating with %d Users ---\n', current_num_users);

    % Initialize user locations for CoMP and Non-CoMP (separate for independence)
    current_user_locations_comp = area_size_m * rand(current_num_users, 2);
    current_user_locations_noncomp = area_size_m * rand(current_num_users, 2);

    % Data Storage for this user count simulation
    all_scheduled_sinrs_linear_comp_rzf_current = [];
    all_scheduled_sinrs_linear_noncomp_rzf_current = [];

    %% Time slot loop for CoMP RZF (for current_num_users)
    fprintf('  Starting CoMP RZF simulation for %d users (%d time slots)...\n', current_num_users, num_time_slots_to_simulate);

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

        % Channel Generation for current time slot
        channel_h_ts = cell(current_num_users, num_BS);
        instantaneous_channel_gain_ts = zeros(current_num_users, num_BS);
        for i = 1:current_num_users
            for j = 1:num_BS
                user_pos = [current_user_locations_comp(i, :), 0];
                bs_pos = bs_locations(j, :);
                d_ij = norm(user_pos - bs_pos);
                if d_ij < 1e-3, d_ij = 1e-3; end

                if j <= 3 % MBS
                    alpha_NLoS_dB = alpha_NLoS_std_dev_dB * randn(1, 1);
                    path_loss_dB_val = 40 * log10(d_ij / 1000) + (-18 * log10(mbs_height_m) + 21 * log10(fc_MHz) + 80 + alpha_NLoS_dB);
                else % LAP-BS
                    elevation_angle_rad = atan2(bs_locations(j, 3), norm(current_user_locations_comp(i, :) - bs_locations(j, 1:2)));
                    path_loss_free_space_dB = 20 * log10((4 * pi * (fc_MHz * 1e6) * d_ij) / c);
                    path_loss_additional_dB = (sigma_LAP - 2) * (1 - elevation_angle_rad / (pi / 2))^2;
                    path_loss_dB_val = path_loss_free_space_dB + path_loss_additional_dB;
                end
                channel_gain_linear_ts_ij = 10^(-path_loss_dB_val / 10);

                if j <= 3 % MBS (Rayleigh fading)
                    channel_h_ts{i, j} = sqrt(channel_gain_linear_ts_ij) * (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
                else % LAP-BS (Rician fading)
                    h_tilde = ones(num_antennas_per_BS, 1);
                    h_bar = (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
                    channel_h_ts{i, j} = sqrt(channel_gain_linear_ts_ij) * (sqrt(rician_K / (1 + rician_K)) * h_tilde + sqrt(1 / (1 + rician_K)) * h_bar);
                end
                instantaneous_channel_gain_ts(i, j) = norm(channel_h_ts{i, j})^2;
            end
        end

        %% User Scheduling (CoMP)
        candidate_sets_ts = cell(1, num_BS);
        for j = 1:num_BS
            selected_users_for_bs_j = [];
            available_users = 1:current_num_users;
            current_channel_gains_for_bs_j = zeros(1, current_num_users);
            for u_idx = 1:current_num_users
                current_channel_gains_for_bs_j(u_idx) = norm(channel_h_ts{u_idx, j})^2;
            end
            [~, initial_user_idx_in_available] = max(current_channel_gains_for_bs_j(available_users));
            initial_user = available_users(initial_user_idx_in_available);
            selected_users_for_bs_j = [selected_users_for_bs_j, initial_user];
            available_users(initial_user_idx_in_available) = [];

            while length(selected_users_for_bs_j) < sus_candidate_set_size && ~isempty(available_users)
                min_max_correlation = inf;
                best_candidate_user = -1;
                best_candidate_user_idx_in_available = -1;
                for u_candidate_idx = 1:length(available_users)
                    u_candidate = available_users(u_candidate_idx);
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
                        best_candidate_user_idx_in_available = u_candidate_idx;
                    end
                end
                if best_candidate_user ~= -1
                    selected_users_for_bs_j = [selected_users_for_bs_j, best_candidate_user];
                    available_users(best_candidate_user_idx_in_available) = [];
                else
                    break;
                end
            end
            candidate_sets_ts{j} = selected_users_for_bs_j;
        end

        all_candidate_users_ts_comp = unique(cell2mat(candidate_sets_ts));
        overlap_users_ts_comp = [];
        user_candidate_counts_ts_comp = zeros(1, current_num_users);
        for i = 1:current_num_users
            for j = 1:num_BS
                if ismember(i, candidate_sets_ts{j})
                    user_candidate_counts_ts_comp(i) = user_candidate_counts_ts_comp(i) + 1;
                end
            end
            if user_candidate_counts_ts_comp(i) > 1
                overlap_users_ts_comp = [overlap_users_ts_comp, i];
            end
        end

        resolved_assigned_users_ts_comp = zeros(1, current_num_users);
        resolved_user_bs_assignment_ts_comp = zeros(1, current_num_users);
        for u = overlap_users_ts_comp
            bs_considering_u = [];
            for j = 1:num_BS
                if ismember(u, candidate_sets_ts{j})
                    bs_considering_u = [bs_considering_u, j];
                end
            end
            max_gain = -inf;
            best_bs_for_u = -1;
            for j = bs_considering_u
                if norm(channel_h_ts{u, j})^2 > max_gain
                    max_gain = norm(channel_h_ts{u, j})^2;
                    best_bs_for_u = j;
                end
            end
            resolved_assigned_users_ts_comp(u) = 1;
            resolved_user_bs_assignment_ts_comp(u) = best_bs_for_u;
        end
        non_overlap_users_ts_comp = setdiff(all_candidate_users_ts_comp, overlap_users_ts_comp);
        for u = non_overlap_users_ts_comp
            for j = 1:num_BS
                if ismember(u, candidate_sets_ts{j})
                    resolved_assigned_users_ts_comp(u) = 1;
                    resolved_user_bs_assignment_ts_comp(u) = j;
                    break;
                end
            end
        end

        initially_assigned_users_ts_comp = find(resolved_assigned_users_ts_comp == 1);
        if length(initially_assigned_users_ts_comp) > max_simultaneous_users_comp
            user_gain_from_assigned_bs_ts_comp = zeros(1, length(initially_assigned_users_ts_comp));
            assigned_user_indices_for_gain_sort_ts_comp = initially_assigned_users_ts_comp;
            for u_idx = 1:length(assigned_user_indices_for_gain_sort_ts_comp)
                u = assigned_user_indices_for_gain_sort_ts_comp(u_idx);
                assigned_bs = resolved_user_bs_assignment_ts_comp(u);
                user_gain_from_assigned_bs_ts_comp(u_idx) = norm(channel_h_ts{u, assigned_bs})^2;
            end
            [~, sorted_indices_in_assigned_ts_comp] = sort(user_gain_from_assigned_bs_ts_comp, 'descend');
            selected_user_indices_in_initially_assigned_ts_comp = sorted_indices_in_assigned_ts_comp(1:max_simultaneous_users_comp);
            S_ts_comp = initially_assigned_users_ts_comp(selected_user_indices_in_initially_assigned_ts_comp);
        else
            S_ts_comp = initially_assigned_users_ts_comp;
        end
        num_scheduled_users_ts_comp = length(S_ts_comp);

        %% RZF Beamforming and SINR Calculation (CoMP)
        if num_scheduled_users_ts_comp > 0
            H_S_ts_comp = zeros(num_scheduled_users_ts_comp, total_antennas);
            for k = 1:num_scheduled_users_ts_comp
                user_idx = S_ts_comp(k);
                row_k = [];
                for j = 1:num_BS
                    row_k = [row_k, channel_h_ts{user_idx, j}'];
                end
                H_S_ts_comp(k, :) = row_k;
            end

            [~, sinrs_comp_rzf_ts] = calculate_rzf_sinr_comp_helper(H_S_ts_comp, P_j_watt, num_antennas_per_BS, noise_power_watt, c_init, c_max, SINR_threshold_linear);
            if ~isempty(sinrs_comp_rzf_ts)
                all_scheduled_sinrs_linear_comp_rzf_current = [all_scheduled_sinrs_linear_comp_rzf_current, sinrs_comp_rzf_ts];
            end
        end
    end
    fprintf('  CoMP RZF simulation complete for %d users.\n', current_num_users);

    %% Time slot loop for Non-CoMP RZF (for current_num_users)
    fprintf('  Starting Non-CoMP RZF simulation for %d users (%d time slots)...\n', current_num_users, num_time_slots_to_simulate);

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

        % Channel Generation for current time slot
        channel_h_ts = cell(current_num_users, num_BS);
        instantaneous_channel_gain_ts = zeros(current_num_users, num_BS);
        for i = 1:current_num_users
            for j = 1:num_BS
                user_pos = [current_user_locations_noncomp(i, :), 0];
                bs_pos = bs_locations(j, :);
                d_ij = norm(user_pos - bs_pos);
                if d_ij < 1e-3, d_ij = 1e-3; end

                if j <= 3 % MBS
                    alpha_NLoS_dB = alpha_NLoS_std_dev_dB * randn(1, 1);
                    path_loss_dB_val = 40 * log10(d_ij / 1000) + (-18 * log10(mbs_height_m) + 21 * log10(fc_MHz) + 80 + alpha_NLoS_dB);
                else % LAP-BS
                    elevation_angle_rad = atan2(bs_locations(j, 3), norm(current_user_locations_noncomp(i, :) - bs_locations(j, 1:2)));
                    path_loss_free_space_dB = 20 * log10((4 * pi * (fc_MHz * 1e6) * d_ij) / c);
                    path_loss_additional_dB = (sigma_LAP - 2) * (1 - elevation_angle_rad / (pi / 2))^2;
                    path_loss_dB_val = path_loss_free_space_dB + path_loss_additional_dB;
                end
                channel_gain_linear_ts_ij = 10^(-path_loss_dB_val / 10);

                if j <= 3 % MBS (Rayleigh fading)
                    channel_h_ts{i, j} = sqrt(channel_gain_linear_ts_ij) * (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
                else % LAP-BS (Rician fading)
                    h_tilde = ones(num_antennas_per_BS, 1);
                    h_bar = (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
                    channel_h_ts{i, j} = sqrt(channel_gain_linear_ts_ij) * (sqrt(rician_K / (1 + rician_K)) * h_tilde + sqrt(1 / (1 + rician_K)) * h_bar);
                end
                instantaneous_channel_gain_ts(i, j) = norm(channel_h_ts{i, j})^2;
            end
        end

        %% User Scheduling (Non-CoMP)
        scheduled_users_ts_noncomp = [];
        scheduled_users_per_bs_ts_noncomp = cell(1, num_BS);
        assigned_bs_for_scheduled_users_noncomp_ts = zeros(1, current_num_users);

        candidate_sets_noncomp_ts = cell(1, num_BS);
        for j = 1:num_BS
            selected_users_for_bs_j_sus = [];
            available_users_for_sus = 1:current_num_users;
            current_channel_gains_for_bs_j_sus = zeros(1, current_num_users);
            for u_idx = 1:current_num_users
                current_channel_gains_for_bs_j_sus(u_idx) = norm(channel_h_ts{u_idx, j})^2;
            end
            [~, initial_user_idx_in_available_sus] = max(current_channel_gains_for_bs_j_sus(available_users_for_sus));
            initial_user_sus = available_users_for_sus(initial_user_idx_in_available_sus);
            selected_users_for_bs_j_sus = [selected_users_for_bs_j_sus, initial_user_sus];
            available_users_for_sus(initial_user_idx_in_available_sus) = [];

            while length(selected_users_for_bs_j_sus) < sus_candidate_set_size && ~isempty(available_users_for_sus)
                min_max_correlation_sus = inf;
                best_candidate_user_sus = -1;
                best_candidate_user_idx_in_available_sus = -1;
                for u_candidate_idx_sus = 1:length(available_users_for_sus)
                    u_candidate_sus = available_users_for_sus(u_candidate_idx_sus);
                    max_correlation_sus = 0;
                    for u_selected_sus = selected_users_for_bs_j_sus
                        norm_u_candidate_sus = norm(channel_h_ts{u_candidate_sus, j});
                        norm_u_selected_sus = norm(channel_h_ts{u_selected_sus, j});
                        if norm_u_candidate_sus > 0 && norm_u_selected_sus > 0
                            correlation_sus = abs(channel_h_ts{u_candidate_sus, j}' * channel_h_ts{u_selected_sus, j}) / (norm_u_candidate_sus * norm_u_selected_sus);
                            max_correlation_sus = max(max_correlation_sus, correlation_sus);
                        else
                            max_correlation_sus = 1;
                        end
                    end
                    if max_correlation_sus < min_max_correlation_sus
                        min_max_correlation_sus = max_correlation_sus;
                        best_candidate_user_sus = u_candidate_sus;
                        best_candidate_user_idx_in_available_sus = u_candidate_idx_sus;
                    end
                end
                if best_candidate_user_sus ~= -1
                    selected_users_for_bs_j_sus = [selected_users_for_bs_j_sus, best_candidate_user_sus];
                    available_users_for_sus(best_candidate_user_idx_in_available_sus) = [];
                else
                    break;
                end
            end
            candidate_sets_noncomp_ts{j} = selected_users_for_bs_j_sus;

            candidate_users_for_bs_j = selected_users_for_bs_j_sus;
            num_candidate_users_j = length(candidate_users_for_bs_j);
            if num_candidate_users_j > 0
                gain_to_bs_j_candidates = zeros(1, num_candidate_users_j);
                for u_idx_in_candidate = 1:num_candidate_users_j
                    u = candidate_users_for_bs_j(u_idx_in_candidate);
                    gain_to_bs_j_candidates(u_idx_in_candidate) = instantaneous_channel_gain_ts(u, j);
                end
                [~, sorted_indices_in_candidate] = sort(gain_to_bs_j_candidates, 'descend');
                num_to_schedule_j = min(num_candidate_users_j, max_simultaneous_users_per_bs_noncomp);
                scheduled_users_j_ts = candidate_users_for_bs_j(sorted_indices_in_candidate(1:num_to_schedule_j));
                scheduled_users_ts_noncomp = [scheduled_users_ts_noncomp, scheduled_users_j_ts];
                scheduled_users_per_bs_ts_noncomp{j} = scheduled_users_j_ts;
                assigned_bs_for_scheduled_users_noncomp_ts(scheduled_users_j_ts) = j;
            end
        end
        num_scheduled_users_ts_noncomp = length(scheduled_users_ts_noncomp);

        %% RZF Beamforming and SINR Calculation (Non-CoMP)
        current_user_sinrs_ts_noncomp = [];
        if num_scheduled_users_ts_noncomp > 0
            W_noncomp_rzf_ts = cell(1, num_BS);
            for j = 1:num_BS
                users_in_Sj_noncomp = scheduled_users_per_bs_ts_noncomp{j};
                num_users_in_Sj_noncomp = length(users_in_Sj_noncomp);
                if num_users_in_Sj_noncomp > 0
                    H_j_noncomp_ts = zeros(num_users_in_Sj_noncomp, num_antennas_per_BS);
                    for k = 1:num_users_in_Sj_noncomp
                        user_idx = users_in_Sj_noncomp(k);
                        H_j_noncomp_ts(k, :) = channel_h_ts{user_idx, j}';
                    end
                    [W_j_noncomp_rzf, ~] = calculate_rzf_sinr_noncomp_helper(H_j_noncomp_ts, P_j_watt, noise_power_watt, c_init, c_max, SINR_threshold_linear);
                    W_noncomp_rzf_ts{j} = W_j_noncomp_rzf;
                else
                    W_noncomp_rzf_ts{j} = zeros(num_antennas_per_BS, 0);
                end
            end

            current_user_sinrs_ts_noncomp_temp = zeros(1, num_scheduled_users_ts_noncomp);
            for k = 1:num_scheduled_users_ts_noncomp
                user_idx = scheduled_users_ts_noncomp(k);
                assigned_bs_idx = assigned_bs_for_scheduled_users_noncomp_ts(user_idx);
                if assigned_bs_idx > 0
                    users_at_assigned_bs = scheduled_users_per_bs_ts_noncomp{assigned_bs_idx};
                    user_col_idx_in_W_j = find(users_at_assigned_bs == user_idx);
                    if ~isempty(user_col_idx_in_W_j) && ~isempty(W_noncomp_rzf_ts{assigned_bs_idx})
                        w_k_noncomp = W_noncomp_rzf_ts{assigned_bs_idx}(:, user_col_idx_in_W_j);
                        h_k_assigned_bs = channel_h_ts{user_idx, assigned_bs_idx};
                        desired_signal_power = abs(h_k_assigned_bs' * w_k_noncomp)^2;
                        interference_power = 0;
                        for m = 1:num_BS
                            if m ~= assigned_bs_idx
                                users_at_interfering_bs = scheduled_users_per_bs_ts_noncomp{m};
                                W_m_noncomp_rzf = W_noncomp_rzf_ts{m};
                                if ~isempty(users_at_interfering_bs) && ~isempty(W_m_noncomp_rzf) && all(users_at_interfering_bs <= current_num_users) && all(users_at_interfering_bs > 0)
                                    h_k_interfering_bs = channel_h_ts{user_idx, m};
                                    for l_idx = 1:length(users_at_interfering_bs)
                                        w_l_interfering = W_m_noncomp_rzf(:, l_idx);
                                        interference_power = interference_power + abs(h_k_interfering_bs' * w_l_interfering)^2;
                                    end
                                end
                            end
                        end
                        denominator = interference_power + noise_power_watt + eps;
                        current_user_sinrs_ts_noncomp_temp(k) = desired_signal_power / denominator;
                    else
                        current_user_sinrs_ts_noncomp_temp(k) = NaN;
                    end
                else
                    current_user_sinrs_ts_noncomp_temp(k) = NaN;
                end
            end
            if ~any(isnan(current_user_sinrs_ts_noncomp_temp)) && ~any(isinf(current_user_sinrs_ts_noncomp_temp))
                current_user_sinrs_ts_noncomp = current_user_sinrs_ts_noncomp_temp;
            end
        end

        %% Store Results for Non-CoMP RZF
        if ~isempty(current_user_sinrs_ts_noncomp)
            all_scheduled_sinrs_linear_noncomp_rzf_current = [all_scheduled_sinrs_linear_noncomp_rzf_current, current_user_sinrs_ts_noncomp];
        end
    end
    fprintf('  Non-CoMP RZF simulation complete for %d users.\n', current_num_users);

    %% Calculate Median SINR for the current number of users
    % CoMP RZF
    if ~isempty(all_scheduled_sinrs_linear_comp_rzf_current)
        valid_sinrs_comp = all_scheduled_sinrs_linear_comp_rzf_current(~isnan(all_scheduled_sinrs_linear_comp_rzf_current) & ~isinf(all_scheduled_sinrs_linear_comp_rzf_current) & all_scheduled_sinrs_linear_comp_rzf_current > 0);
        if ~isempty(valid_sinrs_comp)
            median_sinr_comp_rzf_vs_users(user_count_idx) = median(10*log10(valid_sinrs_comp));
        else
            median_sinr_comp_rzf_vs_users(user_count_idx) = NaN; % Store NaN if no valid SINRs
        end
    else
        median_sinr_comp_rzf_vs_users(user_count_idx) = NaN; % Store NaN if no data
    end

    % Non-CoMP RZF
    if ~isempty(all_scheduled_sinrs_linear_noncomp_rzf_current)
        valid_sinrs_noncomp = all_scheduled_sinrs_linear_noncomp_rzf_current(~isnan(all_scheduled_sinrs_linear_noncomp_rzf_current) & ~isinf(all_scheduled_sinrs_linear_noncomp_rzf_current) & all_scheduled_sinrs_linear_noncomp_rzf_current > 0);
        if ~isempty(valid_sinrs_noncomp)
            median_sinr_noncomp_rzf_vs_users(user_count_idx) = median(10*log10(valid_sinrs_noncomp));
        else
            median_sinr_noncomp_rzf_vs_users(user_count_idx) = NaN; % Store NaN if no valid SINRs
        end
    else
        median_sinr_noncomp_rzf_vs_users(user_count_idx) = NaN; % Store NaN if no data
    end
end % End of loop over different numbers of users

fprintf('\n--- Simulation across different numbers of users Complete ---\n');

%% Plot Median SINR vs. Number of Users
fprintf('Generating Median SINR vs. Number of Users Plot...\n');

figure;
hold on;

% Plot CoMP RZF data
plot(num_users_range, median_sinr_comp_rzf_vs_users, 'b-o', 'DisplayName', 'CoMP (RZF)');

% Plot Non-CoMP RZF data
plot(num_users_range, median_sinr_noncomp_rzf_vs_users, 'r-x', 'DisplayName', 'Non-CoMP (RZF)');

xlabel('Number of Users');
ylabel('Median SINR (dB)');
title('Median SINR vs. Number of Users');
legend;
grid on;
hold off;

fprintf('\nPlotting Complete.\n');

%% Helper Function to Calculate RZF Beams and SINR (CoMP)
function [W_normalized, user_sinrs_linear] = calculate_rzf_sinr_comp_helper(H_S, P_j_watt, num_antennas_per_BS, noise_power_watt, c_init, c_max, SINR_threshold_linear)
    [num_scheduled_users, total_antennas] = size(H_S);
    num_BS = total_antennas / num_antennas_per_BS;
    if num_scheduled_users == 0
        W_normalized = []; user_sinrs_linear = []; return;
    end
    current_c = c_init;
    W_normalized_temp = [];
    current_user_sinrs = [];

    while current_c <= c_max
        avg_power_per_user = (num_BS * P_j_watt) / num_scheduled_users;
        delta = current_c * noise_power_watt / avg_power_per_user;
        try
            matrix_to_invert = (H_S * H_S' + delta * eye(num_scheduled_users));
            if size(matrix_to_invert, 1) ~= size(matrix_to_invert, 2) || rank(matrix_to_invert) < size(matrix_to_invert, 1)
                W_unnormalized_temp = H_S' * pinv(matrix_to_invert);
            else
                W_unnormalized_temp = H_S' * inv(matrix_to_invert);
            end
        catch
            W_unnormalized_temp = H_S' * pinv(H_S * H_S' + delta * eye(num_scheduled_users));
        end

        if any(isnan(W_unnormalized_temp(:))) || any(isinf(W_unnormalized_temp(:)))
            W_normalized = []; user_sinrs_linear = []; return;
        end

        W_normalized_temp = zeros(total_antennas, num_scheduled_users);
        valid_normalization = true;
        for j = 1:num_BS
            antenna_indices = ((j-1)*num_antennas_per_BS + 1) : (j*num_antennas_per_BS);
            W_j_unnormalized = W_unnormalized_temp(antenna_indices, :);
            power_at_bs_j = sum(sum(abs(W_j_unnormalized).^2));
            if power_at_bs_j > P_j_watt
                scaling_factor = sqrt(P_j_watt / power_at_bs_j);
                if isnan(scaling_factor) || isinf(scaling_factor)
                    valid_normalization = false; break;
                end
                W_normalized_temp(antenna_indices, :) = W_j_unnormalized * scaling_factor;
            else
                W_normalized_temp(antenna_indices, :) = W_j_unnormalized;
            end
        end

        if ~valid_normalization || any(isnan(W_normalized_temp(:))) || any(isinf(W_normalized_temp(:)))
            W_normalized = []; user_sinrs_linear = []; return;
        end

        current_user_sinrs = zeros(1, num_scheduled_users);
        for k = 1:num_scheduled_users
            w_k = W_normalized_temp(:, k);
            h_k = H_S(k, :)';
            desired_signal_power = abs(h_k' * w_k)^2;
            interference_power = 0;
            for l = 1:num_scheduled_users
                if l ~= k
                    w_l = W_normalized_temp(:, l);
                    interference_power = interference_power + abs(h_k' * w_l)^2;
                end
            end
            denominator = interference_power + noise_power_watt + eps;
            current_user_sinrs(k) = desired_signal_power / denominator;
        end

        min_sinr_achieved = min(current_user_sinrs);
        if any(isnan(current_user_sinrs)) || any(isinf(current_user_sinrs))
            W_normalized = []; user_sinrs_linear = []; return;
        end

        if min_sinr_achieved >= SINR_threshold_linear
            W_normalized = W_normalized_temp;
            user_sinrs_linear = current_user_sinrs;
            return; % RZF converged, return results
        end
        current_c = current_c + 0.1;
    end

    % If loop finishes without convergence, return results from the last iteration if valid
    if ~isempty(current_user_sinrs) && ~any(isnan(current_user_sinrs)) && ~any(isinf(current_user_sinrs))
        W_normalized = W_normalized_temp;
        user_sinrs_linear = current_user_sinrs;
    else
        W_normalized = [];
        user_sinrs_linear = [];
    end
end

%% Helper Function to Calculate RZF Beams and SINR (Non-CoMP - Per BS)
function [W_normalized, user_sinrs_linear] = calculate_rzf_sinr_noncomp_helper(H_j, P_j_watt, noise_power_watt, c_init, c_max, SINR_threshold_linear)
    [num_users_in_Sj, num_antennas_per_BS] = size(H_j);
    if num_users_in_Sj == 0
        W_normalized = zeros(num_antennas_per_BS, 0); user_sinrs_linear = []; return;
    end
    current_c = c_init;
    W_normalized_temp = [];
    current_user_sinrs_local = [];

    while current_c <= c_max
        avg_power_per_user_j = P_j_watt / num_users_in_Sj;
        delta_j = current_c * noise_power_watt / avg_power_per_user_j;
        try
            matrix_to_invert_j = (H_j * H_j' + delta_j * eye(num_users_in_Sj));
            if size(matrix_to_invert_j, 1) ~= size(matrix_to_invert_j, 2) || rank(matrix_to_invert_j) < size(matrix_to_invert_j, 1)
                W_unnormalized_temp = H_j' * pinv(matrix_to_invert_j);
            else
                W_unnormalized_temp = H_j' * inv(matrix_to_invert_j);
            end
        catch
            W_unnormalized_temp = H_j' * pinv(H_j * H_j' + delta_j * eye(num_users_in_Sj));
        end

        if any(isnan(W_unnormalized_temp(:))) || any(isinf(W_unnormalized_temp(:)))
            W_normalized = zeros(num_antennas_per_BS, 0); user_sinrs_linear = []; return;
        end

        W_normalized_temp = W_unnormalized_temp;
        power_at_bs_j = sum(sum(abs(W_normalized_temp).^2));
        if power_at_bs_j > P_j_watt
            scaling_factor = sqrt(P_j_watt / power_at_bs_j);
            if isnan(scaling_factor) || isinf(scaling_factor)
                W_normalized = zeros(num_antennas_per_BS, 0); user_sinrs_linear = []; return;
            end
            W_normalized_temp = W_normalized_temp * scaling_factor;
        end

        if any(isnan(W_normalized_temp(:))) || any(isinf(W_normalized_temp(:)))
            W_normalized = zeros(num_antennas_per_BS, 0); user_sinrs_linear = []; return;
        end

        current_user_sinrs_local = zeros(1, num_users_in_Sj);
        for k = 1:num_users_in_Sj
            w_k_j = W_normalized_temp(:, k);
            h_k_j = H_j(k, :)';
            desired_signal_power_local = abs(h_k_j' * w_k_j)^2;
            intra_cell_interference_local = 0;
            for l = 1:num_users_in_Sj
                if l ~= k
                    w_l_j = W_normalized_temp(:, l);
                    intra_cell_interference_local = intra_cell_interference_local + abs(h_k_j' * w_l_j)^2;
                end
            end
            denominator_local = intra_cell_interference_local + noise_power_watt + eps;
            current_user_sinrs_local(k) = desired_signal_power_local / denominator_local;
        end

        min_sinr_achieved_local = min(current_user_sinrs_local);
        if any(isnan(current_user_sinrs_local)) || any(isinf(current_user_sinrs_local))
            W_normalized = zeros(num_antennas_per_BS, 0); user_sinrs_linear = []; return;
        end

        if min_sinr_achieved_local >= SINR_threshold_linear
            W_normalized = W_normalized_temp;
            user_sinrs_linear = current_user_sinrs_local;
            return; % RZF converged for this BS
        end
        current_c = current_c + 0.1;
    end

    % If loop finishes without convergence, return results from the last iteration if valid
    if ~isempty(W_normalized_temp) && ~any(isnan(W_normalized_temp(:))) && ~any(isinf(W_normalized_temp(:)))
        W_normalized = W_normalized_temp;
        user_sinrs_linear = current_user_sinrs_local;
    else
        W_normalized = zeros(num_antennas_per_BS, 0);
        user_sinrs_linear = [];
    end
end