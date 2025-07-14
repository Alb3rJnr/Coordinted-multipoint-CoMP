% MATLAB Simulation to compare CoMP RZF, ZF, and MRT SINR variation with increasing number of users

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

%% Simulation Parameters (Varying Number of Users and Fixed Time Slots)
num_users_range = [10, 20, 30, 40, 50]; % Range of user numbers to simulate
num_time_slots_to_simulate = 1000; % Fixed number of time slots for each user count

%% Data Storage for SINR vs. Number of Users
median_sinr_comp_rzf_vs_users = zeros(1, length(num_users_range));
median_sinr_comp_zf_vs_users = zeros(1, length(num_users_range));
median_sinr_comp_mrt_vs_users = zeros(1, length(num_users_range));

%% Simulation Loop (Iterate over different numbers of users)
fprintf('Starting Simulation across different numbers of users for CoMP Beamforming Comparison...\n');

for user_count_idx = 1:length(num_users_range)
    current_num_users = num_users_range(user_count_idx);
    fprintf('\n--- Simulating with %d Users ---\n', current_num_users);

    % Initialize user locations for CoMP
    current_user_locations = area_size_m * rand(current_num_users, 2);

    % Determine mobile and stationary users (60% mobile, 40% stationary)
    num_mobile_users = round(0.6 * current_num_users);
    mobile_user_indices = randperm(current_num_users, num_mobile_users);
    stationary_user_indices = setdiff(1:current_num_users, mobile_user_indices);

    % Data Storage for this user count simulation
    all_scheduled_sinrs_linear_comp_rzf_current = [];
    all_scheduled_sinrs_linear_comp_zf_current = [];
    all_scheduled_sinrs_linear_comp_mrt_current = [];

    %% Time slot loop for CoMP (for current_num_users)
    fprintf('  Starting CoMP simulation for %d users (%d time slots)...\n', current_num_users, num_time_slots_to_simulate);

    for time_slot = 1:num_time_slots_to_simulate
        % Update User Locations (Random Walk Model for mobile users only)
        theta_i = 2 * pi * rand(current_num_users, 1); % Random direction in [0, 2*pi]
        v_i = v_min + (v_max - v_min) * rand(current_num_users, 1); % Random speed in [v_min, v_max]
        for i = mobile_user_indices % Update only mobile users
            x_new = current_user_locations(i, 1) + (v_i(i) / v_max) * D_max * cos(theta_i(i));
            y_new = current_user_locations(i, 2) + (v_i(i) / v_max) * D_max * sin(theta_i(i));
            % Boundary reflection
            if x_new < 0
                x_new = -x_new;
                theta_i(i) = pi - theta_i(i);
            elseif x_new > area_size_m
                x_new = 2 * area_size_m - x_new;
                theta_i(i) = pi - theta_i(i);
            end
            if y_new < 0
                y_new = -y_new;
                theta_i(i) = -theta_i(i);
            elseif y_new > area_size_m
                y_new = 2 * area_size_m - y_new;
                theta_i(i) = -theta_i(i);
            end
            current_user_locations(i, :) = [x_new, y_new];
        end
        current_user_locations(:, 1) = max(0, min(area_size_m, current_user_locations(:, 1)));
        current_user_locations(:, 2) = max(0, min(area_size_m, current_user_locations(:, 2)));

        % Channel Generation for current time slot
        channel_h_ts = cell(current_num_users, num_BS);
        instantaneous_channel_gain_ts = zeros(current_num_users, num_BS);
        for i = 1:current_num_users
            for j = 1:num_BS
                user_pos = [current_user_locations(i, :), 0];
                bs_pos = bs_locations(j, :);
                d_ij = norm(user_pos - bs_pos);
                if d_ij < 1e-3, d_ij = 1e-3; end

                if j <= 3 % MBS
                    alpha_NLoS_dB = alpha_NLoS_std_dev_dB * randn(1, 1);
                    path_loss_dB_val = 40 * log10(d_ij / 1000) + (-18 * log10(mbs_height_m) + 21 * log10(fc_MHz) + 80 + alpha_NLoS_dB);
                else % LAP-BS
                    elevation_angle_rad = atan2(bs_locations(j, 3), norm(current_user_locations(i, :) - bs_locations(j, 1:2)));
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

        %% User Scheduling (CoMP: SUS, Max Gain Overlap Resolution, Top Gain Selection)
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

        %% Beamforming and SINR Calculation for RZF, ZF, and MRT (CoMP Scenario)
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

            % --- RZF Beamforming (CoMP) ---
            [~, sinrs_comp_rzf_ts] = calculate_rzf_sinr_comp_helper(H_S_ts_comp, P_j_watt, num_antennas_per_BS, noise_power_watt, c_init, c_max, SINR_threshold_linear);
            if ~isempty(sinrs_comp_rzf_ts)
                all_scheduled_sinrs_linear_comp_rzf_current = [all_scheduled_sinrs_linear_comp_rzf_current, sinrs_comp_rzf_ts];
            end

            % --- ZF Beamforming (CoMP) ---
            if num_scheduled_users_ts_comp > 0 && num_scheduled_users_ts_comp <= total_antennas
                [~, sinrs_comp_zf_ts] = calculate_zf_sinr_comp_helper(H_S_ts_comp, P_j_watt, num_antennas_per_BS, noise_power_watt);
                if ~isempty(sinrs_comp_zf_ts)
                    all_scheduled_sinrs_linear_comp_zf_current = [all_scheduled_sinrs_linear_comp_zf_current, sinrs_comp_zf_ts];
                end
            end

            % --- MRT Beamforming (CoMP) ---
            [~, sinrs_comp_mrt_ts] = calculate_mrt_sinr_comp_helper(H_S_ts_comp, P_j_watt, num_antennas_per_BS, noise_power_watt);
            if ~isempty(sinrs_comp_mrt_ts)
                all_scheduled_sinrs_linear_comp_mrt_current = [all_scheduled_sinrs_linear_comp_mrt_current, sinrs_comp_mrt_ts];
            end
        end
    end
    fprintf('  CoMP simulation complete for %d users.\n', current_num_users);

    %% Calculate Median SINR for the current number of users
    % CoMP RZF
    if ~isempty(all_scheduled_sinrs_linear_comp_rzf_current)
        valid_sinrs = all_scheduled_sinrs_linear_comp_rzf_current(~isnan(all_scheduled_sinrs_linear_comp_rzf_current) & ~isinf(all_scheduled_sinrs_linear_comp_rzf_current) & all_scheduled_sinrs_linear_comp_rzf_current > 0);
        if ~isempty(valid_sinrs)
            median_sinr_comp_rzf_vs_users(user_count_idx) = median(10*log10(valid_sinrs));
        else
            median_sinr_comp_rzf_vs_users(user_count_idx) = NaN;
        end
    else
        median_sinr_comp_rzf_vs_users(user_count_idx) = NaN;
    end

    % CoMP ZF
    if ~isempty(all_scheduled_sinrs_linear_comp_zf_current)
        valid_sinrs = all_scheduled_sinrs_linear_comp_zf_current(~isnan(all_scheduled_sinrs_linear_comp_zf_current) & ~isinf(all_scheduled_sinrs_linear_comp_zf_current) & all_scheduled_sinrs_linear_comp_zf_current > 0);
        if ~isempty(valid_sinrs)
            median_sinr_comp_zf_vs_users(user_count_idx) = median(10*log10(valid_sinrs));
        else
            median_sinr_comp_zf_vs_users(user_count_idx) = NaN;
        end
    else
        median_sinr_comp_zf_vs_users(user_count_idx) = NaN;
    end

    % CoMP MRT
    if ~isempty(all_scheduled_sinrs_linear_comp_mrt_current)
        valid_sinrs = all_scheduled_sinrs_linear_comp_mrt_current(~isnan(all_scheduled_sinrs_linear_comp_mrt_current) & ~isinf(all_scheduled_sinrs_linear_comp_mrt_current) & all_scheduled_sinrs_linear_comp_mrt_current > 0);
        if ~isempty(valid_sinrs)
            median_sinr_comp_mrt_vs_users(user_count_idx) = median(10*log10(valid_sinrs));
        else
            median_sinr_comp_mrt_vs_users(user_count_idx) = NaN;
        end
    else
        median_sinr_comp_mrt_vs_users(user_count_idx) = NaN;
    end
end % End of loop over different numbers of users

fprintf('\n--- Simulation across different numbers of users Complete ---\n');

%% Plot Median SINR vs. Number of Users (CoMP Beamforming Comparison)
fprintf('Generating Median SINR vs. Number of Users Plot (CoMP Beamforming Comparison)...\n');

figure;
hold on;

% Plot CoMP RZF data
plot(num_users_range, median_sinr_comp_rzf_vs_users, 'b-o', 'DisplayName', 'CoMP (RZF)');

% Plot CoMP ZF data
plot(num_users_range, median_sinr_comp_zf_vs_users, 'g-x', 'DisplayName', 'CoMP (ZF)');

% Plot CoMP MRT data
plot(num_users_range, median_sinr_comp_mrt_vs_users, 'm-*', 'DisplayName', 'CoMP (MRT)');

xlabel('Number of Users');
ylabel('Median SINR (dB)');
title('Median SINR vs. Number of Users with 60% Mobile Users (CoMP Beamforming)');
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

        if any(isnan(current_user_sinrs)) || any(isinf(current_user_sinrs))
            current_c = current_c + 0.1; continue;
        end

        min_sinr_achieved = min(current_user_sinrs);
        if min_sinr_achieved >= SINR_threshold_linear
            W_normalized = W_normalized_temp;
            user_sinrs_linear = current_user_sinrs;
            return;
        end
        current_c = current_c + 0.1;
    end

    if ~isempty(current_user_sinrs) && ~any(isnan(current_user_sinrs)) && ~any(isinf(current_user_sinrs))
        W_normalized = W_normalized_temp;
        user_sinrs_linear = current_user_sinrs;
    else
        W_normalized = [];
        user_sinrs_linear = [];
    end
end

%% Helper Function to Calculate ZF Beams and SINR (CoMP)
function [W_normalized, user_sinrs_linear] = calculate_zf_sinr_comp_helper(H_S, P_j_watt, num_antennas_per_BS, noise_power_watt)
    [num_scheduled_users, total_antennas] = size(H_S);
    num_BS = total_antennas / num_antennas_per_BS;
    if num_scheduled_users == 0 || num_scheduled_users > total_antennas
        W_normalized = []; user_sinrs_linear = []; return;
    end

    try
        matrix_to_invert = H_S * H_S';
        if rcond(matrix_to_invert) < eps
            matrix_to_invert = matrix_to_invert + 1e-12 * eye(num_scheduled_users);
        end
        W_zf_unnormalized = H_S' * pinv(matrix_to_invert);
    catch
        W_zf_unnormalized = H_S' * pinv(H_S * H_S');
    end

    if any(isnan(W_zf_unnormalized(:))) || any(isinf(W_zf_unnormalized(:)))
        W_normalized = []; user_sinrs_linear = []; return;
    end

    W_normalized = zeros(total_antennas, num_scheduled_users);
    valid_normalization = true;
    for j = 1:num_BS
        antenna_indices = ((j-1)*num_antennas_per_BS + 1) : (j*num_antennas_per_BS);
        W_j_unnormalized = W_zf_unnormalized(antenna_indices, :);
        power_at_bs_j = sum(sum(abs(W_j_unnormalized).^2));
        if power_at_bs_j > P_j_watt
            scaling_factor = sqrt(P_j_watt / power_at_bs_j);
            if isnan(scaling_factor) || isinf(scaling_factor)
                valid_normalization = false; break;
            end
            W_normalized(antenna_indices, :) = W_j_unnormalized * scaling_factor;
        else
            W_normalized(antenna_indices, :) = W_j_unnormalized;
        end
    end

    if ~valid_normalization || any(isnan(W_normalized(:))) || any(isinf(W_normalized(:)))
        W_normalized = []; user_sinrs_linear = []; return;
    end

    user_sinrs_linear = zeros(1, num_scheduled_users);
    for k = 1:num_scheduled_users
        w_k = W_normalized(:, k);
        h_k = H_S(k, :)';
        desired_signal_power = abs(h_k' * w_k)^2;
        interference_power = 0;
        for l = 1:num_scheduled_users
            if l ~= k
                w_l = W_normalized(:, l);
                interference_power = interference_power + abs(h_k' * w_l)^2;
            end
        end
        denominator = interference_power + noise_power_watt + eps;
        user_sinrs_linear(k) = desired_signal_power / denominator;
    end

    valid_sinr_indices = ~isnan(user_sinrs_linear) & ~isinf(user_sinrs_linear);
    user_sinrs_linear = user_sinrs_linear(valid_sinr_indices);
    W_normalized = W_normalized(:, valid_sinr_indices);
end

%% Helper Function to Calculate MRT Beams and SINR (CoMP)
function [W_normalized, user_sinrs_linear] = calculate_mrt_sinr_comp_helper(H_S, P_j_watt, num_antennas_per_BS, noise_power_watt)
    [num_scheduled_users, total_antennas] = size(H_S);
    num_BS = total_antennas / num_antennas_per_BS;
    if num_scheduled_users == 0
        W_normalized = []; user_sinrs_linear = []; return;
    end

    W_mrt_unnormalized = H_S';
    if any(isnan(W_mrt_unnormalized(:))) || any(isinf(W_mrt_unnormalized(:)))
        W_normalized = []; user_sinrs_linear = []; return;
    end

    W_normalized = zeros(total_antennas, num_scheduled_users);
    valid_normalization = true;
    for j = 1:num_BS
        antenna_indices = ((j-1)*num_antennas_per_BS + 1) : (j*num_antennas_per_BS);
        W_j_unnormalized = W_mrt_unnormalized(antenna_indices, :);
        power_at_bs_j = sum(sum(abs(W_j_unnormalized).^2));
        if power_at_bs_j > P_j_watt
            scaling_factor = sqrt(P_j_watt / power_at_bs_j);
            if isnan(scaling_factor) || isinf(scaling_factor)
                valid_normalization = false; break;
            end
            W_normalized(antenna_indices, :) = W_j_unnormalized * scaling_factor;
        else
            W_normalized(antenna_indices, :) = W_j_unnormalized;
        end
    end

    if ~valid_normalization || any(isnan(W_normalized(:))) || any(isinf(W_normalized(:)))
        W_normalized = []; user_sinrs_linear = []; return;
    end

    user_sinrs_linear = zeros(1, num_scheduled_users);
    for k = 1:num_scheduled_users
        w_k = W_normalized(:, k);
        h_k = H_S(k, :)';
        desired_signal_power = abs(h_k' * w_k)^2;
        interference_power = 0;
        for l = 1:num_scheduled_users
            if l ~= k
                w_l = W_normalized(:, l);
                interference_power = interference_power + abs(h_k' * w_l)^2;
            end
        end
        denominator = interference_power + noise_power_watt + eps;
        user_sinrs_linear(k) = desired_signal_power / denominator;
    end

    valid_sinr_indices = ~isnan(user_sinrs_linear) & ~isinf(user_sinrs_linear);
    user_sinrs_linear = user_sinrs_linear(valid_sinr_indices);
    W_normalized = W_normalized(:, valid_sinr_indices);
end
