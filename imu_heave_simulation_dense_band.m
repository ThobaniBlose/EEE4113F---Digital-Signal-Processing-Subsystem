%% IMU Heave Simulation - Step 1
% Purpose:
% Create the simulation time vector and the true heave motion:
% displacement eta(t), velocity v(t), and acceleration a(t)

clear; clc; close all;

%% Simulation settings
fs = 100;                  % Sampling frequency [Hz]
Ts = 1/fs;                 % Sampling period [s]
T_total = 30*60;           % Total simulation time [s] = 30 minutes
t = (0:Ts:T_total-Ts)';    % Time vector [s]

%% True heave motion parameters - dense random sea state
% Modelling choice:
% build a denser sea state using many sinusoidal components
% with random phases and a smooth amplitude envelope

N = 30;                           % number of components
f_min_sea = 0.08;                 % minimum wave frequency [Hz]
f_max_sea = 0.40;                 % maximum wave frequency [Hz]
f_vec = linspace(f_min_sea, f_max_sea, N)';   % frequency vector [Hz]

% Smooth amplitude envelope (Gaussian-shaped around a peak frequency)
f_peak_shape = 0.16;              % centre of wave-energy envelope [Hz]
sigma_f = 0.06;                   % spread of envelope [Hz]
A_shape = exp(-0.5*((f_vec - f_peak_shape)/sigma_f).^2);

% Set a target significant wave height for the simulated sea
Hm0_target = 0.30;                % [m]
m0_target = (Hm0_target/4)^2;     % target spectral zeroth moment [m^2]

% Scale amplitudes so that 0.5*sum(A_i^2) = m0_target
A_vec = sqrt(2*m0_target) * A_shape / norm(A_shape);

rng(10);                          % repeatable random phases
phi_vec = 2*pi*rand(N,1);         % random phases [rad]

%% True heave motion - dense random sea state
eta = zeros(size(t));             % displacement [m]
v   = zeros(size(t));             % velocity [m/s]
a   = zeros(size(t));             % acceleration [m/s^2]

for i = 1:N
    eta = eta + A_vec(i)*sin(2*pi*f_vec(i)*t + phi_vec(i));
    v   = v   + 2*pi*f_vec(i)*A_vec(i)*cos(2*pi*f_vec(i)*t + phi_vec(i));
    a   = a   - (2*pi*f_vec(i))^2*A_vec(i)*sin(2*pi*f_vec(i)*t + phi_vec(i));
end
%% Plot true motion
t_zoom_end = 20;                          % zoom window [s]
idx_zoom = t <= t_zoom_end;

figure;
plot(t, eta, 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('\eta(t) [m]');
title('True Heave Displacement - Full Duration');

figure;
plot(t(idx_zoom), eta(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('\eta(t) [m]');
title('True Heave Displacement - First 20 s');

figure;
plot(t, v, 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('v(t) [m/s]');
title('True Heave Velocity - Full Duration');

figure;
plot(t(idx_zoom), v(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('v(t) [m/s]');
title('True Heave Velocity - First 20 s');

figure;
plot(t, a, 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('a(t) [m/s^2]');
title('True Heave Acceleration - Full Duration');

figure;
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('a(t) [m/s^2]');
title('True Heave Acceleration - First 20 s');

%% Ideal IMU outputs for heave-only motion
% Heave-only assumptions:
% - no roll, pitch, or yaw
% - body frame stays aligned with navigation frame
% - accelerometer measures specific force
% - gyroscope measures angular rate

g = 9.81;    % gravitational acceleration [m/s^2]

% Ideal accelerometer output [m/s^2]
% ax = 0, ay = 0, az = g + true heave acceleration
accel_ideal_x = zeros(size(t));
accel_ideal_y = zeros(size(t));
accel_ideal_z = g + a;

% Ideal gyroscope output [rad/s]
gyro_ideal_x = zeros(size(t));
gyro_ideal_y = zeros(size(t));
gyro_ideal_z = zeros(size(t));

%% Plot ideal IMU outputs - first 20 s
figure;
plot(t(idx_zoom), accel_ideal_x(idx_zoom), 'LineWidth', 1.2); hold on;
plot(t(idx_zoom), accel_ideal_y(idx_zoom), 'LineWidth', 1.2);
plot(t(idx_zoom), accel_ideal_z(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('Ideal Accelerometer Output - First 20 s');
legend('a_x', 'a_y', 'a_z', 'Location', 'best');

figure;
plot(t(idx_zoom), gyro_ideal_x(idx_zoom), 'LineWidth', 1.2); hold on;
plot(t(idx_zoom), gyro_ideal_y(idx_zoom), 'LineWidth', 1.2);
plot(t(idx_zoom), gyro_ideal_z(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Angular Rate [rad/s]');
title('Ideal Gyroscope Output - First 20 s');
legend('\omega_x', '\omega_y', '\omega_z', 'Location', 'best');

%% Add constant sensor bias
% Modelling choice for now:
% - constant bias on each accelerometer axis
% - constant bias on each gyroscope axis

% Accelerometer bias [m/s^2]
b_ax = 0.02;
b_ay = -0.01;
b_az = 0.03;

% Gyroscope bias [rad/s]
b_gx = 0.005;
b_gy = -0.003;
b_gz = 0.002;

% Biased accelerometer signals
accel_bias_x = accel_ideal_x + b_ax;
accel_bias_y = accel_ideal_y + b_ay;
accel_bias_z = accel_ideal_z + b_az;

% Biased gyroscope signals
gyro_bias_x = gyro_ideal_x + b_gx;
gyro_bias_y = gyro_ideal_y + b_gy;
gyro_bias_z = gyro_ideal_z + b_gz;

%% Plot biased accelerometer output - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), accel_bias_x(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('a_x [m/s^2]');
title('Biased Accelerometer Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), accel_bias_y(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('a_y [m/s^2]');

subplot(3,1,3);
plot(t(idx_zoom), accel_bias_z(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('a_z [m/s^2]');
xlabel('Time [s]');

%% Plot biased gyroscope output - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), gyro_bias_x(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('\omega_x [rad/s]');
title('Biased Gyroscope Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), gyro_bias_y(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('\omega_y [rad/s]');

subplot(3,1,3);
plot(t(idx_zoom), gyro_bias_z(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('\omega_z [rad/s]');
xlabel('Time [s]');

%% Add white Gaussian noise
% Modelling choice for now:
% - zero-mean
% - white
% - Gaussian
% - independent on each axis

rng(1);   % for repeatability

% Accelerometer noise standard deviation [m/s^2]
sigma_ax = 0.01;
sigma_ay = 0.01;
sigma_az = 0.01;

% Gyroscope noise standard deviation [rad/s]
sigma_gx = 0.002;
sigma_gy = 0.002;
sigma_gz = 0.002;

% Generate accelerometer noise
n_ax = sigma_ax * randn(size(t));
n_ay = sigma_ay * randn(size(t));
n_az = sigma_az * randn(size(t));

% Generate gyroscope noise
n_gx = sigma_gx * randn(size(t));
n_gy = sigma_gy * randn(size(t));
n_gz = sigma_gz * randn(size(t));

% Noisy accelerometer signals
accel_noisy_x = accel_bias_x + n_ax;
accel_noisy_y = accel_bias_y + n_ay;
accel_noisy_z = accel_bias_z + n_az;

% Noisy gyroscope signals
gyro_noisy_x = gyro_bias_x + n_gx;
gyro_noisy_y = gyro_bias_y + n_gy;
gyro_noisy_z = gyro_bias_z + n_gz;

%% Plot noisy accelerometer output - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), accel_noisy_x(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('a_x [m/s^2]');
title('Noisy Accelerometer Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), accel_noisy_y(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('a_y [m/s^2]');

subplot(3,1,3);
plot(t(idx_zoom), accel_noisy_z(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('a_z [m/s^2]');
xlabel('Time [s]');

%% Plot noisy gyroscope output - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), gyro_noisy_x(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('\omega_x [rad/s]');
title('Noisy Gyroscope Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), gyro_noisy_y(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('\omega_y [rad/s]');

subplot(3,1,3);
plot(t(idx_zoom), gyro_noisy_z(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('\omega_z [rad/s]');
xlabel('Time [s]');

%% Apply first-order low-pass filter
% For now, use a simple first-order equivalent LPF.
% These cutoffs are temporary debugging values, not final MPU6050 settings.

fc_accel = 5;      % accelerometer LPF cutoff [Hz]
fc_gyro  = 5;      % gyroscope LPF cutoff [Hz]

alpha_accel = exp(-2*pi*fc_accel*Ts);
alpha_gyro  = exp(-2*pi*fc_gyro*Ts);

% Preallocate filtered signals
accel_filt_x = zeros(size(t));
accel_filt_y = zeros(size(t));
accel_filt_z = zeros(size(t));

gyro_filt_x = zeros(size(t));
gyro_filt_y = zeros(size(t));
gyro_filt_z = zeros(size(t));

% Initial conditions
accel_filt_x(1) = accel_noisy_x(1);
accel_filt_y(1) = accel_noisy_y(1);
accel_filt_z(1) = accel_noisy_z(1);

gyro_filt_x(1) = gyro_noisy_x(1);
gyro_filt_y(1) = gyro_noisy_y(1);
gyro_filt_z(1) = gyro_noisy_z(1);

% Recursive filtering
for k = 2:length(t)
    accel_filt_x(k) = alpha_accel*accel_filt_x(k-1) + (1-alpha_accel)*accel_noisy_x(k);
    accel_filt_y(k) = alpha_accel*accel_filt_y(k-1) + (1-alpha_accel)*accel_noisy_y(k);
    accel_filt_z(k) = alpha_accel*accel_filt_z(k-1) + (1-alpha_accel)*accel_noisy_z(k);

    gyro_filt_x(k) = alpha_gyro*gyro_filt_x(k-1) + (1-alpha_gyro)*gyro_noisy_x(k);
    gyro_filt_y(k) = alpha_gyro*gyro_filt_y(k-1) + (1-alpha_gyro)*gyro_noisy_y(k);
    gyro_filt_z(k) = alpha_gyro*gyro_filt_z(k-1) + (1-alpha_gyro)*gyro_noisy_z(k);
end

%% Plot filtered vs noisy signals - first 20 s
figure;
subplot(2,1,1);
plot(t(idx_zoom), accel_noisy_z(idx_zoom), 'LineWidth', 0.8); hold on;
plot(t(idx_zoom), accel_filt_z(idx_zoom), 'LineWidth', 1.4);
grid on;
ylabel('a_z [m/s^2]');
title('Accelerometer Z-Axis: Noisy vs Filtered');
legend('Noisy', 'Filtered', 'Location', 'best');

subplot(2,1,2);
plot(t(idx_zoom), gyro_noisy_x(idx_zoom), 'LineWidth', 0.8); hold on;
plot(t(idx_zoom), gyro_filt_x(idx_zoom), 'LineWidth', 1.4);
grid on;
xlabel('Time [s]');
ylabel('\omega_x [rad/s]');
title('Gyroscope X-Axis: Noisy vs Filtered');
legend('Noisy', 'Filtered', 'Location', 'best');

%% Plot filtered accelerometer outputs - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), accel_filt_x(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('a_x [m/s^2]');
title('Filtered Accelerometer Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), accel_filt_y(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('a_y [m/s^2]');

subplot(3,1,3);
plot(t(idx_zoom), accel_filt_z(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('a_z [m/s^2]');
xlabel('Time [s]');

%% Plot filtered gyroscope outputs - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), gyro_filt_x(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('\omega_x [rad/s]');
title('Filtered Gyroscope Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), gyro_filt_y(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('\omega_y [rad/s]');

subplot(3,1,3);
plot(t(idx_zoom), gyro_filt_z(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('\omega_z [rad/s]');
xlabel('Time [s]');

%% Apply saturation using chosen full-scale ranges
% Chosen ranges for this simulation:
% Accelerometer: +/-4 g
% Gyroscope: +/-250 deg/s

accel_fs_g = 4;                    % accelerometer full-scale [g]
gyro_fs_dps = 250;                 % gyroscope full-scale [deg/s]

accel_limit = accel_fs_g * g;      % accelerometer saturation limit [m/s^2]
gyro_limit  = deg2rad(gyro_fs_dps);% gyroscope saturation limit [rad/s]

% Saturated accelerometer signals
accel_sat_x = min(max(accel_filt_x, -accel_limit), accel_limit);
accel_sat_y = min(max(accel_filt_y, -accel_limit), accel_limit);
accel_sat_z = min(max(accel_filt_z, -accel_limit), accel_limit);

% Saturated gyroscope signals
gyro_sat_x = min(max(gyro_filt_x, -gyro_limit), gyro_limit);
gyro_sat_y = min(max(gyro_filt_y, -gyro_limit), gyro_limit);
gyro_sat_z = min(max(gyro_filt_z, -gyro_limit), gyro_limit);

%% Plot saturated outputs - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), accel_sat_x(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('a_x [m/s^2]');
title('Saturated Accelerometer Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), accel_sat_y(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('a_y [m/s^2]');

subplot(3,1,3);
plot(t(idx_zoom), accel_sat_z(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('a_z [m/s^2]');
xlabel('Time [s]');

figure;
subplot(3,1,1);
plot(t(idx_zoom), gyro_sat_x(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('\omega_x [rad/s]');
title('Saturated Gyroscope Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), gyro_sat_y(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('\omega_y [rad/s]');

subplot(3,1,3);
plot(t(idx_zoom), gyro_sat_z(idx_zoom), 'LineWidth', 1.1);
grid on;
ylabel('\omega_z [rad/s]');
xlabel('Time [s]');

%% Check whether any clipping occurred
num_clip_accel = sum(abs(accel_filt_x) >= accel_limit) + ...
                 sum(abs(accel_filt_y) >= accel_limit) + ...
                 sum(abs(accel_filt_z) >= accel_limit);

num_clip_gyro = sum(abs(gyro_filt_x) >= gyro_limit) + ...
                sum(abs(gyro_filt_y) >= gyro_limit) + ...
                sum(abs(gyro_filt_z) >= gyro_limit);

fprintf('Total clipped accelerometer samples: %d\n', num_clip_accel);
fprintf('Total clipped gyroscope samples: %d\n', num_clip_gyro);

%% Quantisation using MPU-6050 sensitivity scale factors
% Datasheet-based sensitivity values for the chosen ranges:
% Accelerometer +/-4 g  -> 8192 LSB/g
% Gyroscope    +/-250 dps -> 131 LSB/(deg/s)

accel_sens_lsb_per_g = 8192;   % [LSB/g]
gyro_sens_lsb_per_dps = 131;   % [LSB/(deg/s)]

% Convert saturated accelerometer signals from m/s^2 to g
accel_sat_x_g = accel_sat_x / g;
accel_sat_y_g = accel_sat_y / g;
accel_sat_z_g = accel_sat_z / g;

% Convert saturated gyroscope signals from rad/s to deg/s
gyro_sat_x_dps = rad2deg(gyro_sat_x);
gyro_sat_y_dps = rad2deg(gyro_sat_y);
gyro_sat_z_dps = rad2deg(gyro_sat_z);

% Quantise to integer sensor counts
accel_counts_x = round(accel_sat_x_g * accel_sens_lsb_per_g);
accel_counts_y = round(accel_sat_y_g * accel_sens_lsb_per_g);
accel_counts_z = round(accel_sat_z_g * accel_sens_lsb_per_g);

gyro_counts_x = round(gyro_sat_x_dps * gyro_sens_lsb_per_dps);
gyro_counts_y = round(gyro_sat_y_dps * gyro_sens_lsb_per_dps);
gyro_counts_z = round(gyro_sat_z_dps * gyro_sens_lsb_per_dps);

% Convert quantised counts back to engineering units
accel_quant_x = (accel_counts_x / accel_sens_lsb_per_g) * g;
accel_quant_y = (accel_counts_y / accel_sens_lsb_per_g) * g;
accel_quant_z = (accel_counts_z / accel_sens_lsb_per_g) * g;

gyro_quant_x = deg2rad(gyro_counts_x / gyro_sens_lsb_per_dps);
gyro_quant_y = deg2rad(gyro_counts_y / gyro_sens_lsb_per_dps);
gyro_quant_z = deg2rad(gyro_counts_z / gyro_sens_lsb_per_dps);

%% Quantisation step sizes
accel_lsb_mps2 = g / accel_sens_lsb_per_g;                 % [m/s^2 per LSB]
gyro_lsb_rps   = deg2rad(1 / gyro_sens_lsb_per_dps);       % [rad/s per LSB]

fprintf('Accelerometer quantisation step: %.6e m/s^2 per LSB\n', accel_lsb_mps2);
fprintf('Gyroscope quantisation step: %.6e rad/s per LSB\n', gyro_lsb_rps);

%% Plot quantised vs saturated signals - first 20 s
figure;
subplot(2,1,1);
plot(t(idx_zoom), accel_sat_z(idx_zoom), 'LineWidth', 1.0); hold on;
plot(t(idx_zoom), accel_quant_z(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('a_z [m/s^2]');
title('Accelerometer Z-Axis: Saturated vs Quantised');
legend('Saturated', 'Quantised', 'Location', 'best');

subplot(2,1,2);
plot(t(idx_zoom), gyro_sat_x(idx_zoom), 'LineWidth', 1.0); hold on;
plot(t(idx_zoom), gyro_quant_x(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('\omega_x [rad/s]');
title('Gyroscope X-Axis: Saturated vs Quantised');
legend('Saturated', 'Quantised', 'Location', 'best');

%% Plot quantised accelerometer outputs - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), accel_quant_x(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('a_x [m/s^2]');
title('Quantised Accelerometer Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), accel_quant_y(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('a_y [m/s^2]');

subplot(3,1,3);
plot(t(idx_zoom), accel_quant_z(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('a_z [m/s^2]');
xlabel('Time [s]');

%% Plot quantised gyroscope outputs - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), gyro_quant_x(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('\omega_x [rad/s]');
title('Quantised Gyroscope Output - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), gyro_quant_y(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('\omega_y [rad/s]');

subplot(3,1,3);
plot(t(idx_zoom), gyro_quant_z(idx_zoom), 'LineWidth', 1.0);
grid on;
ylabel('\omega_z [rad/s]');
xlabel('Time [s]');

%% Recover estimated heave acceleration from simulated accelerometer
% Heave-only assumption:
% body frame stays aligned with navigation frame,
% so vertical acceleration estimate is az - g

a_est = accel_quant_z - g;    % estimated heave acceleration [m/s^2]

%% Plot true vs estimated heave acceleration - first 20 s
figure;
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.3); hold on;
plot(t(idx_zoom), a_est(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('True vs Estimated Heave Acceleration - First 20 s');
legend('True a(t)', 'Estimated a_{est}(t)', 'Location', 'best');

%% Plot estimation error - first 20 s
a_err = a_est - a;

figure;
plot(t(idx_zoom), a_err(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Error [m/s^2]');
title('Heave Acceleration Estimation Error - First 20 s');

%% Print basic error statistics
fprintf('Mean acceleration estimation error: %.6e m/s^2\n', mean(a_err));
fprintf('Std acceleration estimation error:  %.6e m/s^2\n', std(a_err));
fprintf('Max abs acceleration error:         %.6e m/s^2\n', max(abs(a_err)));

%% Bias correction of recovered heave acceleration
% For this simulation test, assume the z-axis accelerometer bias is known

a_est_bc = accel_quant_z - g - b_az;   % bias-corrected acceleration estimate

%% Plot true vs bias-corrected estimated acceleration - first 20 s
figure;
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.3); hold on;
plot(t(idx_zoom), a_est_bc(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('True vs Bias-Corrected Heave Acceleration - First 20 s');
legend('True a(t)', 'Bias-corrected a_{est}(t)', 'Location', 'best');

%% Plot bias-corrected estimation error - first 20 s
a_err_bc = a_est_bc - a;

figure;
plot(t(idx_zoom), a_err_bc(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Error [m/s^2]');
title('Bias-Corrected Heave Acceleration Error - First 20 s');

%% Print bias-corrected error statistics
fprintf('Mean bias-corrected acceleration error: %.6e m/s^2\n', mean(a_err_bc));
fprintf('Std bias-corrected acceleration error:  %.6e m/s^2\n', std(a_err_bc));
fprintf('Max abs bias-corrected acceleration error: %.6e m/s^2\n', max(abs(a_err_bc)));

%% PSD of true and recovered heave acceleration using Welch's method
% This is a validation step on acceleration first.
% Later we can move to displacement spectrum / wave spectrum.

% Welch settings (chosen for this simulation)
win_len = 4096;                  % window length [samples]
noverlap = win_len/2;            % 50% overlap
nfft = 8192;                     % FFT length

% Remove any small DC offset before PSD estimation
a_true_psd = detrend(a, 'constant');
a_est_psd  = detrend(a_est_bc, 'constant');

% Welch PSD estimates
[Paa_true, f_psd] = pwelch(a_true_psd, win_len, noverlap, nfft, fs);
[Paa_est,  ~]     = pwelch(a_est_psd,  win_len, noverlap, nfft, fs);

%% Plot PSD comparison
figure;
plot(f_psd, Paa_true, 'LineWidth', 1.4); hold on;
plot(f_psd, Paa_est, 'LineWidth', 1.2);
grid on;
xlabel('Frequency [Hz]');
ylabel('PSD [(m/s^2)^2/Hz]');
title('Acceleration PSD: True vs Bias-Corrected Estimate');
legend('True acceleration PSD', 'Estimated acceleration PSD', 'Location', 'best');
xlim([0 2]);   % focus on low-frequency wave band

%% Find dominant frequency in the low-frequency band
f_min = 0.02;
f_max = 2.0;
idx_band = (f_psd >= f_min) & (f_psd <= f_max);

[~, idx_true_peak_local] = max(Paa_true(idx_band));
[~, idx_est_peak_local]  = max(Paa_est(idx_band));

f_band = f_psd(idx_band);

f_true_peak = f_band(idx_true_peak_local);
f_est_peak  = f_band(idx_est_peak_local);

fprintf('True acceleration PSD peak frequency: %.6f Hz\n', f_true_peak);
fprintf('Estimated acceleration PSD peak frequency: %.6f Hz\n', f_est_peak);

%% Convert acceleration PSD to wave-elevation PSD
% Spectral double integration:
% S_eta(f) = S_a(f) / (2*pi*f)^4
%
% Important:
% avoid division near f = 0 by applying a low-frequency cutoff

f_int_min = 0.05;    % low-frequency cutoff [Hz] for safe division

omega_psd = 2*pi*f_psd;

% True elevation PSD directly from eta(t)
eta_true_psd = detrend(eta, 'constant');
[See_true, f_eta] = pwelch(eta_true_psd, win_len, noverlap, nfft, fs);

% Estimated elevation PSD from bias-corrected acceleration PSD
See_est = zeros(size(Paa_est));

idx_valid = f_psd >= f_int_min;
See_est(idx_valid) = Paa_est(idx_valid) ./ (omega_psd(idx_valid).^4);

%% Plot true vs estimated wave-elevation PSD
figure;
plot(f_eta, See_true, 'LineWidth', 1.4); hold on;
plot(f_psd, See_est, 'LineWidth', 1.2);
grid on;
xlabel('Frequency [Hz]');
ylabel('PSD [m^2/Hz]');
title('Wave Elevation PSD: True vs Estimated from Acceleration');
legend('True elevation PSD', 'Estimated elevation PSD', 'Location', 'best');
xlim([0 1]);

%% Find dominant frequency in the valid wave band
f_wave_min = 0.05;
f_wave_max = 1.0;

idx_wave = (f_eta >= f_wave_min) & (f_eta <= f_wave_max);
idx_wave_est = (f_psd >= f_wave_min) & (f_psd <= f_wave_max);

[~, idx_true_eta_peak_local] = max(See_true(idx_wave));
[~, idx_est_eta_peak_local]  = max(See_est(idx_wave_est));

f_true_eta_band = f_eta(idx_wave);
f_est_eta_band  = f_psd(idx_wave_est);

f_true_eta_peak = f_true_eta_band(idx_true_eta_peak_local);
f_est_eta_peak  = f_est_eta_band(idx_est_eta_peak_local);

fprintf('True elevation PSD peak frequency: %.6f Hz\n', f_true_eta_peak);
fprintf('Estimated elevation PSD peak frequency: %.6f Hz\n', f_est_eta_peak);

%% Compute wave parameters from the elevation PSD
% Use the same valid wave band for both true and estimated spectra

f_param_min = 0.05;
f_param_max = 1.0;

idx_true_band = (f_eta >= f_param_min) & (f_eta <= f_param_max);
idx_est_band  = (f_psd >= f_param_min) & (f_psd <= f_param_max);

f_true_band = f_eta(idx_true_band);
f_est_band  = f_psd(idx_est_band);

See_true_band = See_true(idx_true_band);
See_est_band  = See_est(idx_est_band);

% Spectral moments for true spectrum
m0_true = trapz(f_true_band, See_true_band);
m1_true = trapz(f_true_band, f_true_band .* See_true_band);
m2_true = trapz(f_true_band, (f_true_band.^2) .* See_true_band);

% Spectral moments for estimated spectrum
m0_est = trapz(f_est_band, See_est_band);
m1_est = trapz(f_est_band, f_est_band .* See_est_band);
m2_est = trapz(f_est_band, (f_est_band.^2) .* See_est_band);

% Peak frequency and peak period
[~, idx_pk_true] = max(See_true_band);
[~, idx_pk_est]  = max(See_est_band);

fp_true = f_true_band(idx_pk_true);
fp_est  = f_est_band(idx_pk_est);

Tp_true = 1 / fp_true;
Tp_est  = 1 / fp_est;

% Derived wave parameters
Hm0_true  = 4 * sqrt(m0_true);
Hm0_est   = 4 * sqrt(m0_est);

Tm01_true = m0_true / m1_true;
Tm01_est  = m0_est / m1_est;

Tm02_true = sqrt(m0_true / m2_true);
Tm02_est  = sqrt(m0_est / m2_est);

%% Print wave parameters
fprintf('\n--- True spectrum wave parameters ---\n');
fprintf('m0   = %.6e m^2\n', m0_true);
fprintf('Hm0  = %.6f m\n', Hm0_true);
fprintf('fp   = %.6f Hz\n', fp_true);
fprintf('Tp   = %.6f s\n', Tp_true);
fprintf('Tm01 = %.6f s\n', Tm01_true);
fprintf('Tm02 = %.6f s\n', Tm02_true);

fprintf('\n--- Estimated spectrum wave parameters ---\n');
fprintf('m0   = %.6e m^2\n', m0_est);
fprintf('Hm0  = %.6f m\n', Hm0_est);
fprintf('fp   = %.6f Hz\n', fp_est);
fprintf('Tp   = %.6f s\n', Tp_est);
fprintf('Tm01 = %.6f s\n', Tm01_est);
fprintf('Tm02 = %.6f s\n', Tm02_est);

%% Compare estimated elevation PSD with the true input component frequencies
% This helps verify that the recovered spectrum contains the frequencies
% that were intentionally placed in the simulated sea state.

figure;
plot(f_eta, See_true, 'LineWidth', 1.4); hold on;
plot(f_psd, See_est, 'LineWidth', 1.2);

% Mark the input component frequencies
for i = 1:length(f_vec)
    xline(f_vec(i), '--k', sprintf('f_%d = %.2f Hz', i, f_vec(i)), ...
        'LabelVerticalAlignment', 'middle', ...
        'LabelHorizontalAlignment', 'left');
end

grid on;
xlabel('Frequency [Hz]');
ylabel('PSD [m^2/Hz]');
title('Wave Elevation PSD with Input Component Frequencies');
legend('True elevation PSD', 'Estimated elevation PSD', 'Location', 'best');
xlim([0 0.5]);

%% Inspect the input sea-state components
figure;
stem(f_vec, A_vec, 'filled');
grid on;
xlabel('Frequency [Hz]');
ylabel('Amplitude A_i [m]');
title('Input Multi-Sine Sea-State Components');
xlim([0.05 0.42]);