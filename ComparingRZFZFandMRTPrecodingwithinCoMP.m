% MATLAB Simulation for Coordinated Beamforming (CoMP)
% Comparing RZF, ZF, and MRT Precoding within CoMP

% Clear workspace and command window
clear;
clc;
close all;

%% System Parameters
num_users = 50; % Total number of users (I)
num_BS = 4; % Total number of base stations (J) (1 LAP, 3 MBS)
num_antennas_per_BS = 4; % Number of transmit antennas per BS (Nt)
total_antennas = num_BS * num_antennas_per_BS; % Total antennas in the cluster

% Simulation area (square)
area_size_m = 5000; % 5 km x 5 km

% Coverage Radii in meters
mbs_coverage_radius_m = 1000; % 1 km
lap_coverage_radius_m = 2000; % 2 km

% BS Locations [x, y, z] in meters
center_x = area_size_m / 2; % 2500
center_y = area_size_m / 2; % 2500
mbs_triangle_side_length = 2 * mbs_coverage_radius_m; % 2000 m
mbs_distance_to_center = mbs_triangle_side_length / sqrt(3);

bs_locations = [
    center_x, center_y + mbs_distance_to_center, 30; % MBS 1
    center_x - mbs_triangle_side_length/2, center_y - mbs_distance_to_center/2, 30; % MBS 2
    center_x + mbs_triangle_side_length/2, center_y - mbs_distance_to_center/2, 30; % MBS 3
    center_x, center_y, 2000; % LAP-BS
];

% User Locations [x, y] in meters - Randomly distributed
user_locations = area_size_m * rand(num_users, 2);

% Initialize mobility status (60% mobile, 40% stationary)
num_mobile_users = floor(0.6 * num_users); % 60% of users are mobile
mobile_user_indices = randperm(num_users, num_mobile_users); % Randomly select mobile users
is_mobile_user = zeros(1, num_users); % 1 for mobile, 0 for stationary
is_mobile_user(mobile_user_indices) = 1;

% Random Walk Mobility Parameters
v_min = 1; % Minimum speed (m/s)
v_max = 10; % Maximum speed (m/s)
D_max = 10; % Maximum distance per time slot (meters)

% Channel Model Parameters
fc_MHz = 2545; % Carrier frequency (2 GHz)
c = 3e8; % Speed of light
lambda = c / (fc_MHz * 1e6); % Wavelength
noise_power_dBm = -100; % Noise power
noise_power_watt = 10^(noise_power_dBm / 10) * 1e-3; % Noise power in Watts

% Path Loss Model Parameters
mbs_height_m = 30; % MBS height
alpha_NLoS_std_dev_dB = 6; % Shadowing standard deviation
sigma_LAP = 4; % Elevation-dependent loss for LAP
rician_K_dB = 10; % Rician K factor in dB
rician_K = 10^(rician_K_dB / 10); % Rician K factor (linear)

% Power Constraints
P_j_dBm = 46; % Per-BS transmit power (1 Watt)
P_j_watt = 10^(P_j_dBm / 10) * 1e-3; % Per-BS power in Watts
P_total_watt = num_BS * P_j_watt; % Total power

% Scheduling & Beamforming Parameters
sus_candidate_set_size = 8; % Candidate users per BS
SINR_threshold_dB = 0; % Minimum SINR for RZF
SINR_threshold_linear = 10^(SINR_threshold_dB / 10);
c_init = 1; % Initial RZF regularization
c_max = 100; % Maximum RZF regularization
max_simultaneous_users_comp = total_antennas; %

%% Simulation Parameters
num_time_slots_to_simulate = 1000; % Number of time slots

%% Data Storage for Multi-Time Slot Simulation
all_scheduled_sinrs_linear_comp_rzf = [];
all_scheduled_sinrs_linear_comp_zf = [];
all_scheduled_sinrs_linear_comp_mrt = [];

%% Network Topology Visualization (Initial, static)
figure;
hold on;
scatter(bs_locations(1:3, 1), bs_locations(1:3, 2), 100, 's', 'filled', 'DisplayName', 'MBS');
scatter(bs_locations(4, 1), bs_locations(4, 2), 150, '^', 'filled', 'DisplayName', 'LAP-BS');
scatter(user_locations(:, 1), user_locations(:, 2), 30, 'o', 'filled', 'DisplayName', 'Users');
theta = 0:0.01:2*pi;
for j = 1:3
    x_circle = bs_locations(j, 1) + mbs_coverage_radius_m * cos(theta);
    y_circle = bs_locations(j, 2) + mbs_coverage_radius_m * sin(theta);
    plot(x_circle, y_circle, 'b--', 'DisplayName', sprintf('MBS %d Coverage', j));
end
x_circle = bs_locations(4, 1) + lap_coverage_radius_m * cos(theta);
y_circle = bs_locations(4, 2) + lap_coverage_radius_m * sin(theta);
plot(x_circle, y_circle, 'g--', 'DisplayName', 'LAP-BS Coverage');
xlabel('X-coordinate (m)');
ylabel('Y-coordinate (m)');
title('Network Topology: BS, User Locations, and Coverage Areas');
legend;
grid on;
axis equal;
xlim([0 area_size_m]);
ylim([0 area_size_m]);
hold off;

%% Multi-Time Slot Simulation Loop
fprintf('Starting Multi-Time Slot Simulation (%d slots) for CoMP scenarios with Random Walk mobility and binary scheduling...\n', num_time_slots_to_simulate);
for time_slot = 1:num_time_slots_to_simulate
    fprintf('\n--- Time Slot %d ---\n', time_slot);
    %% Update User Locations (Random Walk Model for mobile users only)
    % Generate random directions and speeds
    theta_i = 2 * pi * rand(num_users, 1); % Random direction in [0, 2*pi]
    v_i = v_min + (v_max - v_min) * rand(num_users, 1); % Random speed in [v_min, v_max]
    % Update locations
    for i = 1:num_users
        if is_mobile_user(i) % Update only mobile users
            x_new = user_locations(i, 1) + (v_i(i) / v_max) * D_max * cos(theta_i(i));
            y_new = user_locations(i, 2) + (v_i(i) / v_max) * D_max * sin(theta_i(i));
            % Boundary reflection
            if x_new < 0
                x_new = -x_new; % Reflect off x=0
                theta_i(i) = pi - theta_i(i); % Reverse x-component of direction
            elseif x_new > area_size_m
                x_new = 2 * area_size_m - x_new; % Reflect off x=area_size_m
                theta_i(i) = pi - theta_i(i);
            end
            if y_new < 0
                y_new = -y_new; % Reflect off y=0
                theta_i(i) = -theta_i(i); % Reverse y-component of direction
            elseif y_new > area_size_m
                y_new = 2 * area_size_m - y_new; % Reflect off y=area_size_m
                theta_i(i) = -theta_i(i);
            end
            user_locations(i, :) = [x_new, y_new];
        end
    end
    % Ensure users stay within bounds (additional check)
    user_locations(:, 1) = max(0, min(area_size_m, user_locations(:, 1)));
    user_locations(:, 2) = max(0, min(area_size_m, user_locations(:, 2)));

    %% Channel Generation (Large-scale and small-scale fading)
    % Large-scale fading (path loss and shadowing)
    channel_gain_linear = zeros(num_users, num_BS);
    for i = 1:num_users
        for j = 1:num_BS
            user_pos = [user_locations(i, :), 0];
            bs_pos = bs_locations(j, :);
            d_ij = norm(user_pos - bs_pos);
            if d_ij < 1e-3
                d_ij = 1e-3;
            end
            if j <= 3 % MBS
                alpha_NLoS_dB = alpha_NLoS_std_dev_dB * randn(1, 1);
                path_loss_dB_val = 40 * log10(d_ij / 1000) + (-18 * log10(mbs_height_m) + 21 * log10(fc_MHz) + 80 + alpha_NLoS_dB);
            else % LAP-BS
                elevation_angle_rad = atan2(bs_locations(j, 3), norm(user_locations(i, :) - bs_locations(j, 1:2)));
                path_loss_free_space_dB = 20 * log10((4 * pi * (fc_MHz * 1e6) * d_ij) / c);
                path_loss_additional_dB = (sigma_LAP - 2) * (1 - elevation_angle_rad / (pi / 2))^2;
                path_loss_dB_val = path_loss_free_space_dB + path_loss_additional_dB;
            end
            channel_gain_linear(i, j) = 10^(-path_loss_dB_val / 10);
        end
    end
    % Small-scale fading
    channel_h_ts = cell(num_users, num_BS);
    for i = 1:num_users
        for j = 1:num_BS
            if j <= 3 % MBS (Rayleigh)
                channel_h_ts{i, j} = sqrt(channel_gain_linear(i, j)) * (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
            else % LAP-BS (Rician)
                h_tilde = ones(num_antennas_per_BS, 1);
                h_bar = (randn(num_antennas_per_BS, 1) + 1i * randn(num_antennas_per_BS, 1)) / sqrt(2);
                channel_h_ts{i, j} = sqrt(channel_gain_linear(i, j)) * (sqrt(rician_K / (1 + rician_K)) * h_tilde + sqrt(1 / (1 + rician_K)) * h_bar);
            end
        end
    end
    % Construct full channel matrix H for all users
    H_ts_comp = zeros(num_users, total_antennas);
    for i = 1:num_users
        row_i = [];
        for j = 1:num_BS
            row_i = [row_i, channel_h_ts{i, j}'];
        end
        H_ts_comp(i, :) = row_i;
    end

    %% User Scheduling (CoMP SUS)
    candidate_sets_ts = cell(1, num_BS);
    for j = 1:num_BS
        selected_users_for_bs_j = [];
        available_users = 1:num_users;
        current_channel_gains_for_bs_j = zeros(1, num_users);
        for u_idx = 1:num_users
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
    all_candidate_users_ts = unique(cell2mat(candidate_sets_ts));
    overlap_users_ts = [];
    user_candidate_counts_ts = zeros(1, num_users);
    for i = 1:num_users
        for j = 1:num_BS
            if ismember(i, candidate_sets_ts{j})
                user_candidate_counts_ts(i) = user_candidate_counts_ts(i) + 1;
            end
        end
        if user_candidate_counts_ts(i) > 1
            overlap_users_ts = [overlap_users_ts, i];
        end
    end
    resolved_assigned_users_ts = zeros(1, num_users);
    resolved_user_bs_assignment_ts = zeros(1, num_users);
    for u = overlap_users_ts
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
        resolved_assigned_users_ts(u) = 1;
        resolved_user_bs_assignment_ts(u) = best_bs_for_u;
    end
    non_overlap_users_ts = setdiff(all_candidate_users_ts, overlap_users_ts);
    for u = non_overlap_users_ts
        for j = 1:num_BS
            if ismember(u, candidate_sets_ts{j})
                resolved_assigned_users_ts(u) = 1;
                resolved_user_bs_assignment_ts(u) = j;
                break;
            end
        end
    end
    initially_assigned_users_ts = find(resolved_assigned_users_ts == 1);
    if length(initially_assigned_users_ts) > max_simultaneous_users_comp
        user_gain_from_assigned_bs_ts = zeros(1, length(initially_assigned_users_ts));
        assigned_user_indices_for_gain_sort_ts = initially_assigned_users_ts;
        for u_idx = 1:length(assigned_user_indices_for_gain_sort_ts)
            u = assigned_user_indices_for_gain_sort_ts(u_idx);
            assigned_bs = resolved_user_bs_assignment_ts(u);
            user_gain_from_assigned_bs_ts(u_idx) = norm(channel_h_ts{u, assigned_bs})^2;
        end
        [~, sorted_indices_in_assigned_ts] = sort(user_gain_from_assigned_bs_ts, 'descend');
        selected_user_indices_in_initially_assigned_ts = sorted_indices_in_assigned_ts(1:max_simultaneous_users_comp);
        scheduled_user_indices_ts_comp = initially_assigned_users_ts(selected_user_indices_in_initially_assigned_ts);
    else
        scheduled_user_indices_ts_comp = initially_assigned_users_ts;
    end
    S_ts_comp = scheduled_user_indices_ts_comp;
    num_scheduled_users_ts_comp = length(S_ts_comp);
    final_user_bs_assignment_ts_comp = resolved_user_bs_assignment_ts(S_ts_comp);
    % Define binary scheduling variable
    s_k = zeros(1, num_users);
    s_k(S_ts_comp) = 1; % 1 for scheduled users, 0 otherwise

    %% Beamforming and SINR Calculation
    if num_scheduled_users_ts_comp > 0
        % RZF
        [W_comp_rzf_ts, sinrs_comp_rzf_ts] = calculate_rzf_sinr(H_ts_comp, s_k, S_ts_comp, P_j_watt, num_antennas_per_BS, noise_power_watt, c_init, c_max, SINR_threshold_linear);
        if ~isempty(sinrs_comp_rzf_ts)
            all_scheduled_sinrs_linear_comp_rzf = [all_scheduled_sinrs_linear_comp_rzf, sinrs_comp_rzf_ts];
        end
        % ZF
        if num_scheduled_users_ts_comp <= total_antennas
            [W_comp_zf_ts, sinrs_comp_zf_ts] = calculate_zf_sinr(H_ts_comp, s_k, S_ts_comp, P_j_watt, num_antennas_per_BS, noise_power_watt);
            if ~isempty(sinrs_comp_zf_ts)
                all_scheduled_sinrs_linear_comp_zf = [all_scheduled_sinrs_linear_comp_zf, sinrs_comp_zf_ts];
            end
        end
        % MRT
        [W_comp_mrt_ts, sinrs_comp_mrt_ts] = calculate_mrt_sinr(H_ts_comp, s_k, S_ts_comp, P_j_watt, num_antennas_per_BS, noise_power_watt);
        if ~isempty(sinrs_comp_mrt_ts)
            all_scheduled_sinrs_linear_comp_mrt = [all_scheduled_sinrs_linear_comp_mrt, sinrs_comp_mrt_ts];
        end
    else
        fprintf('No users scheduled for CoMP in Time Slot %d.\n', time_slot);
    end
end
fprintf('\n--- Multi-Time Slot Simulation Complete ---\n');
fprintf('Total time slots simulated: %d\n', num_time_slots_to_simulate);

%% SINR CCDF Plot
fprintf('Generating SINR CCDF Comparison Plot (CoMP Only)...\n');
figure;
hold on;
if ~isempty(all_scheduled_sinrs_linear_comp_rzf)
    valid_indices = ~isnan(all_scheduled_sinrs_linear_comp_rzf) & ~isinf(all_scheduled_sinrs_linear_comp_rzf) & all_scheduled_sinrs_linear_comp_rzf > 0;
    sinrs_dB = 10*log10(all_scheduled_sinrs_linear_comp_rzf(valid_indices));
    sorted_sinrs_dB = sort(sinrs_dB);
    num_scheduled_instances = length(sorted_sinrs_dB);
    ccdf_values = 1 - (1:num_scheduled_instances) / num_scheduled_instances + (1/num_scheduled_instances);
    semilogy(sorted_sinrs_dB, ccdf_values, 'b-', 'DisplayName', 'CoMP (RZF)');
end
if ~isempty(all_scheduled_sinrs_linear_comp_zf)
    valid_indices = ~isnan(all_scheduled_sinrs_linear_comp_zf) & ~isinf(all_scheduled_sinrs_linear_comp_zf) & all_scheduled_sinrs_linear_comp_zf > 0;
    sinrs_dB = 10*log10(all_scheduled_sinrs_linear_comp_zf(valid_indices));
    sorted_sinrs_dB = sort(sinrs_dB);
    num_scheduled_instances = length(sorted_sinrs_dB);
    ccdf_values = 1 - (1:num_scheduled_instances) / num_scheduled_instances + (1/num_scheduled_instances);
    semilogy(sorted_sinrs_dB, ccdf_values, 'b--', 'DisplayName', 'CoMP (ZF)');
end
if ~isempty(all_scheduled_sinrs_linear_comp_mrt)
    valid_indices = ~isnan(all_scheduled_sinrs_linear_comp_mrt) & ~isinf(all_scheduled_sinrs_linear_comp_mrt) & all_scheduled_sinrs_linear_comp_mrt > 0;
    sinrs_dB = 10*log10(all_scheduled_sinrs_linear_comp_mrt(valid_indices));
    sorted_sinrs_dB = sort(sinrs_dB);
    num_scheduled_instances = length(sorted_sinrs_dB);
    ccdf_values = 1 - (1:num_scheduled_instances) / num_scheduled_instances + (1/num_scheduled_instances);
    semilogy(sorted_sinrs_dB, ccdf_values, 'b:', 'DisplayName', 'CoMP (MRT)');
end
xlabel('SINR (dB)');
ylabel('CCDF P(SINR > x)');
title('SINR CCDF Comparison (CoMP: RZF, ZF, MRT)');
legend;
grid on;
hold off;

fprintf('\nSimulation and Plotting Complete.\n');

%% Helper Function to Calculate RZF Beams and SINR
function [W_normalized, user_sinrs_linear] = calculate_rzf_sinr(H, s_k, scheduled_user_indices, P_j_watt, num_antennas_per_BS, noise_power_watt, c_init, c_max, SINR_threshold_linear)
    num_users = size(H, 1);
    total_antennas = size(H, 2);
    num_BS = total_antennas / num_antennas_per_BS;
    num_scheduled_users = length(scheduled_user_indices);
    if num_scheduled_users == 0
        W_normalized = [];
        user_sinrs_linear = [];
        return;
    end
    % Select channel matrix for scheduled users
    H_S = H(scheduled_user_indices, :);
    current_c = c_init;
    rzf_converged = false;
    W_normalized = [];
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
            break;
        end
        W_normalized_temp = zeros(total_antennas, num_scheduled_users);
        valid_normalization = true;
        for j = 1:num_BS
            antenna_indices = ((j-1)*num_antennas_per_BS + 1):(j*num_antennas_per_BS);
            W_j_unnormalized = W_unnormalized_temp(antenna_indices, :);
            power_at_bs_j = sum(sum(abs(W_j_unnormalized).^2));
            if power_at_bs_j > P_j_watt
                scaling_factor = sqrt(P_j_watt / power_at_bs_j);
                if isnan(scaling_factor) || isinf(scaling_factor)
                    valid_normalization = false;
                    break;
                end
                W_normalized_temp(antenna_indices, :) = W_j_unnormalized * scaling_factor;
            else
                W_normalized_temp(antenna_indices, :) = W_j_unnormalized;
            end
        end
        if ~valid_normalization || any(isnan(W_normalized_temp(:))) || any(isinf(W_normalized_temp(:)))
            break;
        end
        % Compute SINR for all users, using s_k
        user_sinrs_linear = zeros(1, num_users);
        for k = 1:num_users
            if s_k(k) == 0
                user_sinrs_linear(k) = 0; % No SINR if not scheduled
                continue;
            end
            scheduled_idx = find(scheduled_user_indices == k);
            if isempty(scheduled_idx)
                user_sinrs_linear(k) = 0;
                continue;
            end
            w_k = W_normalized_temp(:, scheduled_idx);
            h_k = H(k, :)';
            desired_signal_power = s_k(k) * abs(h_k' * w_k)^2;
            interference_power = 0;
            for l = 1:num_users
                if l ~= k && s_k(l) == 1
                    scheduled_idx_l = find(scheduled_user_indices == l);
                    if ~isempty(scheduled_idx_l)
                        w_l = W_normalized_temp(:, scheduled_idx_l);
                        interference_power = interference_power + s_k(l) * abs(h_k' * w_l)^2;
                    end
                end
            end
            denominator = interference_power + noise_power_watt + eps;
            user_sinrs_linear(k) = desired_signal_power / denominator;
        end
        % Check SINR for scheduled users only
        scheduled_sinrs = user_sinrs_linear(scheduled_user_indices);
        min_sinr_achieved = min(scheduled_sinrs);
        if any(isnan(scheduled_sinrs)) || any(isinf(scheduled_sinrs))
            break;
        end
        if min_sinr_achieved >= SINR_threshold_linear
            rzf_converged = true;
            W_normalized = W_normalized_temp;
            user_sinrs_linear = scheduled_sinrs; % Return only scheduled users' SINRs
            break;
        end
        current_c = current_c + 0.1;
    end
    if ~rzf_converged
        if ~isempty(scheduled_sinrs) && ~any(isnan(scheduled_sinrs)) && ~any(isinf(scheduled_sinrs))
            W_normalized = W_normalized_temp;
            user_sinrs_linear = scheduled_sinrs;
        else
            W_normalized = [];
            user_sinrs_linear = [];
        end
    end
end

%% Helper Function to Calculate ZF Beams and SINR
function [W_normalized, user_sinrs_linear] = calculate_zf_sinr(H, s_k, scheduled_user_indices, P_j_watt, num_antennas_per_BS, noise_power_watt)
    num_users = size(H, 1);
    total_antennas = size(H, 2);
    num_BS = total_antennas / num_antennas_per_BS;
    num_scheduled_users = length(scheduled_user_indices);
    if num_scheduled_users == 0 || num_scheduled_users > total_antennas
        W_normalized = [];
        user_sinrs_linear = [];
        return;
    end
    % Select channel matrix for scheduled users
    H_S = H(scheduled_user_indices, :);
    W_zf_unnormalized = H_S' * pinv(H_S * H_S');
    if any(isnan(W_zf_unnormalized(:))) || any(isinf(W_zf_unnormalized(:)))
        W_normalized = [];
        user_sinrs_linear = [];
        return;
    end
    W_normalized = zeros(total_antennas, num_scheduled_users);
    valid_normalization = true;
    for j = 1:num_BS
        antenna_indices = ((j-1)*num_antennas_per_BS + 1):(j*num_antennas_per_BS);
        W_j_unnormalized = W_zf_unnormalized(antenna_indices, :);
        power_at_bs_j = sum(sum(abs(W_j_unnormalized).^2));
        if power_at_bs_j > P_j_watt
            scaling_factor = sqrt(P_j_watt / power_at_bs_j);
            if isnan(scaling_factor) || isinf(scaling_factor)
                valid_normalization = false;
                break;
            end
            W_normalized(antenna_indices, :) = W_j_unnormalized * scaling_factor;
        else
            W_normalized(antenna_indices, :) = W_j_unnormalized;
        end
    end
    if ~valid_normalization || any(isnan(W_normalized(:))) || any(isinf(W_normalized(:)))
        W_normalized = [];
        user_sinrs_linear = [];
        return;
    end
    % Compute SINR for all users, using s_k
    user_sinrs_linear = zeros(1, num_users);
    for k = 1:num_users
        if s_k(k) == 0
            user_sinrs_linear(k) = 0; % No SINR if not scheduled
            continue;
        end
        scheduled_idx = find(scheduled_user_indices == k);
        if isempty(scheduled_idx)
            user_sinrs_linear(k) = 0;
            continue;
        end
        w_k = W_normalized(:, scheduled_idx);
        h_k = H(k, :)';
        desired_signal_power = s_k(k) * abs(h_k' * w_k)^2;
        interference_power = 0;
        for l = 1:num_users
            if l ~= k && s_k(l) == 1
                scheduled_idx_l = find(scheduled_user_indices == l);
                if ~isempty(scheduled_idx_l)
                    w_l = W_normalized(:, scheduled_idx_l);
                    interference_power = interference_power + s_k(l) * abs(h_k' * w_l)^2;
                end
            end
        end
        denominator = interference_power + noise_power_watt + eps;
        user_sinrs_linear(k) = desired_signal_power / denominator;
    end
    user_sinrs_linear = user_sinrs_linear(scheduled_user_indices); % Return only scheduled users' SINRs
    valid_sinr_indices = ~isnan(user_sinrs_linear) & ~isinf(user_sinrs_linear);
    user_sinrs_linear = user_sinrs_linear(valid_sinr_indices);
    W_normalized = W_normalized(:, valid_sinr_indices);
end

%% Helper Function to Calculate MRT Beams and SINR
function [W_normalized, user_sinrs_linear] = calculate_mrt_sinr(H, s_k, scheduled_user_indices, P_j_watt, num_antennas_per_BS, noise_power_watt)
    num_users = size(H, 1);
    total_antennas = size(H, 2);
    num_BS = total_antennas / num_antennas_per_BS;
    num_scheduled_users = length(scheduled_user_indices);
    if num_scheduled_users == 0
        W_normalized = [];
        user_sinrs_linear = [];
        return;
    end
    % Select channel matrix for scheduled users
    H_S = H(scheduled_user_indices, :);
    W_mrt_unnormalized = H_S';
    if any(isnan(W_mrt_unnormalized(:))) || any(isinf(W_mrt_unnormalized(:)))
        W_normalized = [];
        user_sinrs_linear = [];
        return;
    end
    W_normalized = zeros(total_antennas, num_scheduled_users);
    valid_normalization = true;
    for j = 1:num_BS
        antenna_indices = ((j-1)*num_antennas_per_BS + 1):(j*num_antennas_per_BS);
        W_j_unnormalized = W_mrt_unnormalized(antenna_indices, :);
        power_at_bs_j = sum(sum(abs(W_j_unnormalized).^2));
        if power_at_bs_j > P_j_watt
            scaling_factor = sqrt(P_j_watt / power_at_bs_j);
            if isnan(scaling_factor) || isinf(scaling_factor)
                valid_normalization = false;
                break;
            end
            W_normalized(antenna_indices, :) = W_j_unnormalized * scaling_factor;
        else
            W_normalized(antenna_indices, :) = W_j_unnormalized;
        end
    end
    if ~valid_normalization || any(isnan(W_normalized(:))) || any(isinf(W_normalized(:)))
        W_normalized = [];
        user_sinrs_linear = [];
        return;
    end
    % Compute SINR for all users, using s_k
    user_sinrs_linear = zeros(1, num_users);
    for k = 1:num_users
        if s_k(k) == 0
            user_sinrs_linear(k) = 0; % No SINR if not scheduled
            continue;
        end
        scheduled_idx = find(scheduled_user_indices == k);
        if isempty(scheduled_idx)
            user_sinrs_linear(k) = 0;
            continue;
        end
        w_k = W_normalized(:, scheduled_idx);
        h_k = H(k, :)';
        desired_signal_power = s_k(k) * abs(h_k' * w_k)^2;
        interference_power = 0;
        for l = 1:num_users
            if l ~= k && s_k(l) == 1
                scheduled_idx_l = find(scheduled_user_indices == l);
                if ~isempty(scheduled_idx_l)
                    w_l = W_normalized(:, scheduled_idx_l);
                    interference_power = interference_power + s_k(l) * abs(h_k' * w_l)^2;
                end
            end
        end
        denominator = interference_power + noise_power_watt + eps;
        user_sinrs_linear(k) = desired_signal_power / denominator;
    end
    user_sinrs_linear = user_sinrs_linear(scheduled_user_indices); % Return only scheduled users' SINRs
    valid_sinr_indices = ~isnan(user_sinrs_linear) & ~isinf(user_sinrs_linear);
    user_sinrs_linear = user_sinrs_linear(valid_sinr_indices);
    W_normalized = W_normalized(:, valid_sinr_indices);
end
