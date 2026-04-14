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

%% True buoy attitude motion (small roll and pitch)
% Next realism step:
% add small roll and pitch oscillations to the buoy motion

A_roll_deg  = 3.0;          % roll amplitude [deg]
A_pitch_deg = 2.0;          % pitch amplitude [deg]

f_roll  = 0.14;             % roll frequency [Hz]
f_pitch = 0.22;             % pitch frequency [Hz]

phi_roll  = 0.6;            % roll phase [rad]
phi_pitch = -0.4;           % pitch phase [rad]

% Convert amplitudes to radians
A_roll  = deg2rad(A_roll_deg);
A_pitch = deg2rad(A_pitch_deg);

% True roll and pitch angles
roll_true  = A_roll  * sin(2*pi*f_roll*t  + phi_roll);      % [rad]
pitch_true = A_pitch * sin(2*pi*f_pitch*t + phi_pitch);     % [rad]

% True roll and pitch rates
roll_rate_true  = 2*pi*f_roll*A_roll  * cos(2*pi*f_roll*t  + phi_roll);   % [rad/s]
pitch_rate_true = 2*pi*f_pitch*A_pitch * cos(2*pi*f_pitch*t + phi_pitch);  % [rad/s]

% True yaw remains zero for now
yaw_true = zeros(size(t));                 % [rad]
yaw_rate_true = zeros(size(t));            % [rad/s]

%% Plot true attitude motion - first 20 s
figure;
subplot(2,1,1);
plot(t(idx_zoom), rad2deg(roll_true(idx_zoom)), 'LineWidth', 1.2); hold on;
plot(t(idx_zoom), rad2deg(pitch_true(idx_zoom)), 'LineWidth', 1.2);
grid on;
ylabel('Angle [deg]');
title('True Roll and Pitch - First 20 s');
legend('Roll', 'Pitch', 'Location', 'best');

subplot(2,1,2);
plot(t(idx_zoom), rad2deg(roll_rate_true(idx_zoom)), 'LineWidth', 1.2); hold on;
plot(t(idx_zoom), rad2deg(pitch_rate_true(idx_zoom)), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('Rate [deg/s]');
title('True Roll and Pitch Rates - First 20 s');
legend('Roll rate', 'Pitch rate', 'Location', 'best');

%% Ideal IMU outputs with roll and pitch
% Now include the effect of buoy tilt on the accelerometer.
% Navigation frame z-axis is positive upward.
% Specific force in navigation frame is [0; 0; a(t)+g].
% It is then rotated into the body frame.

g = 9.81;    % gravitational acceleration [m/s^2]

% Preallocate ideal accelerometer outputs
accel_ideal_x = zeros(size(t));
accel_ideal_y = zeros(size(t));
accel_ideal_z = zeros(size(t));

% Preallocate ideal gyroscope outputs
gyro_ideal_x = zeros(size(t));
gyro_ideal_y = zeros(size(t));
gyro_ideal_z = zeros(size(t));

for k = 1:length(t)
    % Euler angles at this time
    phi   = roll_true(k);      % roll [rad]
    theta = pitch_true(k);     % pitch [rad]
    psi   = yaw_true(k);       % yaw [rad]

    % Rotation from body to navigation (ZYX convention)
    Rx = [1 0 0;
          0 cos(phi) -sin(phi);
          0 sin(phi)  cos(phi)];

    Ry = [ cos(theta) 0 sin(theta);
           0          1 0;
          -sin(theta) 0 cos(theta)];

    Rz = [cos(psi) -sin(psi) 0;
          sin(psi)  cos(psi) 0;
          0         0        1];

    C_nb = Rz * Ry * Rx;   % body -> navigation
    C_bn = C_nb.';         % navigation -> body

    % Specific force in navigation frame
    f_n = [0; 0; a(k) + g];

    % Ideal accelerometer output in body frame
    f_b = C_bn * f_n;

    accel_ideal_x(k) = f_b(1);
    accel_ideal_y(k) = f_b(2);
    accel_ideal_z(k) = f_b(3);
end

% Ideal gyroscope output in body rates p, q, r
% Use ZYX Euler-angle-rate to body-rate mapping
phi = roll_true;
theta = pitch_true;
psi_dot = yaw_rate_true;

phi_dot = roll_rate_true;
theta_dot = pitch_rate_true;

gyro_ideal_x = phi_dot - sin(theta).*psi_dot;
gyro_ideal_y = cos(phi).*theta_dot + sin(phi).*cos(theta).*psi_dot;
gyro_ideal_z = -sin(phi).*theta_dot + cos(phi).*cos(theta).*psi_dot;

%% Plot ideal IMU outputs - first 20 s
figure;
subplot(3,1,1);
plot(t(idx_zoom), accel_ideal_x(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('a_x [m/s^2]');
title('Ideal Accelerometer Output with Tilt - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), accel_ideal_y(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('a_y [m/s^2]');

subplot(3,1,3);
plot(t(idx_zoom), accel_ideal_z(idx_zoom), 'LineWidth', 1.2);
grid on;
ylabel('a_z [m/s^2]');
xlabel('Time [s]');

figure;
subplot(3,1,1);
plot(t(idx_zoom), rad2deg(gyro_ideal_x(idx_zoom)), 'LineWidth', 1.2);
grid on;
ylabel('p [deg/s]');
title('Ideal Gyroscope Output with Tilt - First 20 s');

subplot(3,1,2);
plot(t(idx_zoom), rad2deg(gyro_ideal_y(idx_zoom)), 'LineWidth', 1.2);
grid on;
ylabel('q [deg/s]');

subplot(3,1,3);
plot(t(idx_zoom), rad2deg(gyro_ideal_z(idx_zoom)), 'LineWidth', 1.2);
grid on;
ylabel('r [deg/s]');
xlabel('Time [s]');

%% Compare naive and tilt-corrected vertical acceleration (ideal case)
% This diagnostic uses the ideal accelerometer and the true roll/pitch.
% It shows why az - g is no longer valid once the buoy tilts.

a_est_naive_ideal = accel_ideal_z - g;   % old heave-only shortcut
a_est_tilt_ideal  = zeros(size(t));      % corrected using true attitude

for k = 1:length(t)
    phi   = roll_true(k);      % roll [rad]
    theta = pitch_true(k);     % pitch [rad]

    % For yaw = 0, the body-to-navigation rotation is Ry*Rx
    Rx = [1 0 0;
          0 cos(phi) -sin(phi);
          0 sin(phi)  cos(phi)];

    Ry = [ cos(theta) 0 sin(theta);
           0          1 0;
          -sin(theta) 0 cos(theta)];

    C_nb = Ry * Rx;          % body -> navigation

    f_b = [accel_ideal_x(k);
           accel_ideal_y(k);
           accel_ideal_z(k)];

    f_n_est = C_nb * f_b;    % rotate specific force back to navigation frame

    % Vertical acceleration estimate = vertical specific force - g contribution
    a_est_tilt_ideal(k) = f_n_est(3) - g;
end

%% Plot true vs naive vs tilt-corrected acceleration - first 20 s
figure;
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.4); hold on;
plot(t(idx_zoom), a_est_naive_ideal(idx_zoom), 'LineWidth', 1.1);
plot(t(idx_zoom), a_est_tilt_ideal(idx_zoom), 'LineWidth', 1.1);
grid on;
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('True vs Naive vs Tilt-Corrected Vertical Acceleration (Ideal Case)');
legend('True a(t)', 'Naive a_z-g', 'Tilt-corrected', 'Location', 'best');

%% Error comparison - first 20 s
err_naive_ideal = a_est_naive_ideal - a;
err_tilt_ideal  = a_est_tilt_ideal  - a;

figure;
plot(t(idx_zoom), err_naive_ideal(idx_zoom), 'LineWidth', 1.1); hold on;
plot(t(idx_zoom), err_tilt_ideal(idx_zoom), 'LineWidth', 1.1);
grid on;
xlabel('Time [s]');
ylabel('Error [m/s^2]');
title('Naive vs Tilt-Corrected Error (Ideal Case)');
legend('Naive error', 'Tilt-corrected error', 'Location', 'best');

%% Print ideal-case error statistics
fprintf('\n--- Ideal tilted case: vertical acceleration recovery ---\n');
fprintf('Naive mean error:          %.6e m/s^2\n', mean(err_naive_ideal));
fprintf('Naive std error:           %.6e m/s^2\n', std(err_naive_ideal));
fprintf('Naive max abs error:       %.6e m/s^2\n', max(abs(err_naive_ideal)));

fprintf('Tilt-corrected mean error: %.6e m/s^2\n', mean(err_tilt_ideal));
fprintf('Tilt-corrected std error:  %.6e m/s^2\n', std(err_tilt_ideal));
fprintf('Tilt-corrected max abs error: %.6e m/s^2\n', max(abs(err_tilt_ideal)));

%% Add sensor bias and z-axis drift
% Accelerometer bias [m/s^2]
b_ax = 0.02;
b_ay = -0.01;

% z-axis accelerometer bias = constant part + random-walk-like drift
b_az0 = 0.03;                 % nominal bias level [m/s^2]

rng(20);                      % repeatable drift realization

b_az_rw = cumsum(randn(size(t)));      % random-walk-like drift
b_az_rw = b_az_rw - mean(b_az_rw);     % remove offset
b_az_rw = 0.005 * b_az_rw / std(b_az_rw);   % set drift strength

b_az = b_az0 + b_az_rw;

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

%% Recover vertical acceleration from non-ideal tilted IMU
% Compare:
% 1) naive method: az - g
% 2) tilt-corrected method using full 3-axis accel and true attitude
% 3) bias-corrected + tilt-corrected method using true injected biases

a_est_naive = accel_quant_z - g;
a_est_tilt = zeros(size(t));
a_est_tilt_bc = zeros(size(t));

for k = 1:length(t)
    phi   = roll_true(k);
    theta = pitch_true(k);
    psi   = yaw_true(k);

    Rx = [1 0 0;
          0 cos(phi) -sin(phi);
          0 sin(phi)  cos(phi)];

    Ry = [ cos(theta) 0 sin(theta);
           0          1 0;
          -sin(theta) 0 cos(theta)];

    Rz = [cos(psi) -sin(psi) 0;
          sin(psi)  cos(psi) 0;
          0         0        1];

    C_nb = Rz * Ry * Rx;   % body -> navigation

    % Uncorrected measured body specific force
    f_b = [accel_quant_x(k);
           accel_quant_y(k);
           accel_quant_z(k)];

    % Bias-corrected body specific force
    f_b_bc = [accel_quant_x(k) - b_ax;
              accel_quant_y(k) - b_ay;
              accel_quant_z(k) - b_az(k)];

    % Rotate to navigation frame
    f_n_est = C_nb * f_b;
    f_n_est_bc = C_nb * f_b_bc;

    % Recover vertical acceleration
    a_est_tilt(k) = f_n_est(3) - g;
    a_est_tilt_bc(k) = f_n_est_bc(3) - g;
end

%% Plot true vs recovered vertical acceleration - first 20 s
figure;
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.3); hold on;
plot(t(idx_zoom), a_est_naive(idx_zoom), 'LineWidth', 1.0);
plot(t(idx_zoom), a_est_tilt(idx_zoom), 'LineWidth', 1.0);
plot(t(idx_zoom), a_est_tilt_bc(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Acceleration [m/s^2]');
title('True vs Recovered Vertical Acceleration (Non-Ideal Tilted Case)');
legend('True a(t)', 'Naive a_z-g', 'Tilt-corrected', 'Bias-corrected + tilt-corrected', ...
       'Location', 'best');

%% Error comparison
err_naive = a_est_naive - a;
err_tilt  = a_est_tilt - a;
err_tilt_bc = a_est_tilt_bc - a;

figure;
plot(t(idx_zoom), err_naive(idx_zoom), 'LineWidth', 1.0); hold on;
plot(t(idx_zoom), err_tilt(idx_zoom), 'LineWidth', 1.0);
plot(t(idx_zoom), err_tilt_bc(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Error [m/s^2]');
title('Recovery Error Comparison (Non-Ideal Tilted Case)');
legend('Naive error', 'Tilt-corrected error', 'Bias-corrected + tilt-corrected error', ...
       'Location', 'best');

%% Print error statistics
fprintf('\n--- Non-ideal tilted case: vertical acceleration recovery ---\n');
fprintf('Naive mean error:                %.6e m/s^2\n', mean(err_naive));
fprintf('Naive std error:                 %.6e m/s^2\n', std(err_naive));
fprintf('Naive max abs error:             %.6e m/s^2\n', max(abs(err_naive)));

fprintf('Tilt-corrected mean error:       %.6e m/s^2\n', mean(err_tilt));
fprintf('Tilt-corrected std error:        %.6e m/s^2\n', std(err_tilt));
fprintf('Tilt-corrected max abs error:    %.6e m/s^2\n', max(abs(err_tilt)));

fprintf('Bias-corrected tilt mean error:  %.6e m/s^2\n', mean(err_tilt_bc));
fprintf('Bias-corrected tilt std error:   %.6e m/s^2\n', std(err_tilt_bc));
fprintf('Bias-corrected tilt max abs error: %.6e m/s^2\n', max(abs(err_tilt_bc)));

%% PSD of true and recovered vertical acceleration (tilted case)
% Use the best recovered signal from the current stage:
% bias-corrected + tilt-corrected vertical acceleration

win_len = 4096;              % window length [samples]
noverlap = win_len/2;        % 50% overlap
nfft = 8192;                 % FFT length

a_true_psd = detrend(a, 'constant');
a_est_psd  = detrend(a_est_tilt_bc, 'constant');

[Paa_true, f_psd] = pwelch(a_true_psd, win_len, noverlap, nfft, fs);
[Paa_est,  ~]     = pwelch(a_est_psd,  win_len, noverlap, nfft, fs);

%% Plot acceleration PSD comparison
figure;
plot(f_psd, Paa_true, 'LineWidth', 1.4); hold on;
plot(f_psd, Paa_est, 'LineWidth', 1.2);
grid on;
xlabel('Frequency [Hz]');
ylabel('PSD [(m/s^2)^2/Hz]');
title('Acceleration PSD: True vs Bias-Corrected Tilt-Corrected Estimate');
legend('True acceleration PSD', 'Estimated acceleration PSD', 'Location', 'best');
xlim([0 2]);

%% Dominant frequency comparison
f_min = 0.02;
f_max = 2.0;
idx_band = (f_psd >= f_min) & (f_psd <= f_max);

[~, idx_true_peak_local] = max(Paa_true(idx_band));
[~, idx_est_peak_local]  = max(Paa_est(idx_band));

f_band = f_psd(idx_band);

f_true_peak = f_band(idx_true_peak_local);
f_est_peak  = f_band(idx_est_peak_local);

fprintf('\n--- Tilted case acceleration PSD ---\n');
fprintf('True acceleration PSD peak frequency: %.6f Hz\n', f_true_peak);
fprintf('Estimated acceleration PSD peak frequency: %.6f Hz\n', f_est_peak);

%% Convert tilted-case acceleration PSD to wave-elevation PSD
% Use the recovered vertical acceleration PSD from the tilted case.
% Temporary working low-frequency cutoff:
f_int_min = 0.08;    % [Hz]

omega_psd = 2*pi*f_psd;

% True elevation PSD directly from eta(t)
eta_true_psd = detrend(eta, 'constant');
[See_true, f_eta] = pwelch(eta_true_psd, win_len, noverlap, nfft, fs);

% Estimated elevation PSD from bias-corrected + tilt-corrected acceleration PSD
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
title('Wave Elevation PSD: True vs Tilt-Corrected Estimate');
legend('True elevation PSD', 'Estimated elevation PSD', 'Location', 'best');
xlim([0.05 0.4]);
ylim([0 0.07]);

%% Find dominant frequency in the useful wave band
f_wave_min = 0.08;
f_wave_max = 1.0;

idx_wave_true = (f_eta >= f_wave_min) & (f_eta <= f_wave_max);
idx_wave_est  = (f_psd >= f_wave_min) & (f_psd <= f_wave_max);

[~, idx_true_eta_peak_local] = max(See_true(idx_wave_true));
[~, idx_est_eta_peak_local]  = max(See_est(idx_wave_est));

f_true_eta_band = f_eta(idx_wave_true);
f_est_eta_band  = f_psd(idx_wave_est);

f_true_eta_peak = f_true_eta_band(idx_true_eta_peak_local);
f_est_eta_peak  = f_est_eta_band(idx_est_eta_peak_local);

fprintf('\n--- Tilted case elevation PSD ---\n');
fprintf('True elevation PSD peak frequency: %.6f Hz\n', f_true_eta_peak);
fprintf('Estimated elevation PSD peak frequency: %.6f Hz\n', f_est_eta_peak);

%% Compute wave parameters from the tilted-case elevation PSD
% Use the same practical band as the current tilted recovery setup.

f_param_min = 0.08;
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

%% Print tilted-case wave parameters
fprintf('\n--- Tilted case true spectrum wave parameters ---\n');
fprintf('m0   = %.6e m^2\n', m0_true);
fprintf('Hm0  = %.6f m\n', Hm0_true);
fprintf('fp   = %.6f Hz\n', fp_true);
fprintf('Tp   = %.6f s\n', Tp_true);
fprintf('Tm01 = %.6f s\n', Tm01_true);
fprintf('Tm02 = %.6f s\n', Tm02_true);

fprintf('\n--- Tilted case estimated spectrum wave parameters ---\n');
fprintf('m0   = %.6e m^2\n', m0_est);
fprintf('Hm0  = %.6f m\n', Hm0_est);
fprintf('fp   = %.6f Hz\n', fp_est);
fprintf('Tp   = %.6f s\n', Tp_est);
fprintf('Tm01 = %.6f s\n', Tm01_est);
fprintf('Tm02 = %.6f s\n', Tm02_est);

%% Percentage errors
Hm0_err_pct  = 100 * (Hm0_est  - Hm0_true)  / Hm0_true;
fp_err_pct   = 100 * (fp_est   - fp_true)   / fp_true;
Tp_err_pct   = 100 * (Tp_est   - Tp_true)   / Tp_true;
Tm01_err_pct = 100 * (Tm01_est - Tm01_true) / Tm01_true;
Tm02_err_pct = 100 * (Tm02_est - Tm02_true) / Tm02_true;

fprintf('\n--- Tilted case wave-parameter errors ---\n');
fprintf('Hm0 error  = %.3f %%\n', Hm0_err_pct);
fprintf('fp error   = %.3f %%\n', fp_err_pct);
fprintf('Tp error   = %.3f %%\n', Tp_err_pct);
fprintf('Tm01 error = %.3f %%\n', Tm01_err_pct);
fprintf('Tm02 error = %.3f %%\n', Tm02_err_pct);

%% Simulated GNSS low-rate vertical position measurement model
% First GNSS model:
% low-rate position-type measurements derived from the true elevation eta(t)

fs_gnss = 1;                      % GNSS update rate [Hz] - safe starting point
Ts_gnss = 1/fs_gnss;
t_gnss = (0:Ts_gnss:T_total-Ts_gnss)';    % GNSS sample times

% True vertical position sampled at GNSS times
eta_true_gnss = interp1(t, eta, t_gnss, 'linear');

% GNSS position noise
rng(40);                          % repeatable noise realization
sigma_eta_gnss = 2.0;             % [m], based on stated positioning precision

% Simulated GNSS vertical position measurement
eta_gnss_meas = eta_true_gnss + sigma_eta_gnss * randn(size(t_gnss));

%% Plot true elevation and GNSS measurements - first 120 s
figure;
plot(t, eta, 'LineWidth', 1.0); hold on;
plot(t_gnss, eta_gnss_meas, 'o', 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Vertical position [m]');
title('True Elevation and Simulated GNSS Position Measurements');
legend('True \eta(t)', 'GNSS samples', 'Location', 'best');
xlim([0 120]);

%% Zoomed plot - first 20 s
figure;
plot(t(t <= 20), eta(t <= 20), 'LineWidth', 1.2); hold on;
plot(t_gnss(t_gnss <= 20), eta_gnss_meas(t_gnss <= 20), 'o', 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Vertical position [m]');
title('True Elevation and Simulated GNSS Position Measurements - First 20 s');
legend('True \eta(t)', 'GNSS samples', 'Location', 'best');

%% GNSS measurement error statistics
eta_gnss_err = eta_gnss_meas - eta_true_gnss;

fprintf('\n--- Simulated GNSS position measurements ---\n');
fprintf('GNSS update rate: %.2f Hz\n', fs_gnss);
fprintf('GNSS position noise std used: %.3f m\n', sigma_eta_gnss);
fprintf('Mean GNSS position error: %.6e m\n', mean(eta_gnss_err));
fprintf('Std GNSS position error:  %.6e m\n', std(eta_gnss_err));
fprintf('Max abs GNSS position error: %.6e m\n', max(abs(eta_gnss_err)));

%% 1D IMU-GNSS Kalman fusion (vertical position and velocity)
% State:
% x = [eta; v]
%
% Process model:
% eta(k+1) = eta(k) + Ts*v(k) + 0.5*Ts^2*u(k)
% v(k+1)   = v(k) + Ts*u(k)
%
% Measurement:
% z = eta_gnss_meas

A_kf = [1 Ts;
        0 1];

B_kf = [0.5*Ts^2;
        Ts];

C_kf = [1 0];

% Initial state estimate
x_hat = [0; 0];

% Initial covariance
P = diag([1, 1]);

% Process-noise tuning
% This is a filter tuning parameter, not a datasheet spec.
sigma_a_proc = 0.05;   % [m/s^2], initial tuning choice
Q = sigma_a_proc^2 * [Ts^4/4, Ts^3/2;
                      Ts^3/2, Ts^2];

% Measurement-noise covariance from GNSS position model
R = sigma_eta_gnss^2;

% Storage
x_hat_hist = zeros(2, length(t));
eta_kf = zeros(size(t));
v_kf   = zeros(size(t));

% GNSS timing
gnss_step = round(fs / fs_gnss);
gnss_idx = 1;

for k = 1:length(t)
    % IMU-driven propagation
    u_k = a_est_tilt_bc(k);

    x_hat = A_kf * x_hat + B_kf * u_k;
    P = A_kf * P * A_kf' + Q;

    % GNSS update when available
    if mod(k-1, gnss_step) == 0 && gnss_idx <= length(t_gnss)
        z_k = eta_gnss_meas(gnss_idx);

        y_k = z_k - C_kf * x_hat;                 % innovation
        S_k = C_kf * P * C_kf' + R;               % innovation covariance
        K_k = P * C_kf' / S_k;                    % Kalman gain

        x_hat = x_hat + K_k * y_k;
        P = (eye(2) - K_k * C_kf) * P;

        gnss_idx = gnss_idx + 1;
    end

    x_hat_hist(:,k) = x_hat;
    eta_kf(k) = x_hat(1);
    v_kf(k)   = x_hat(2);
end

%% Plot fused vertical position - first 120 s
figure;
plot(t, eta, 'LineWidth', 1.2); hold on;
plot(t, eta_kf, 'LineWidth', 1.1);
plot(t_gnss, eta_gnss_meas, 'o');
grid on;
xlabel('Time [s]');
ylabel('Vertical position [m]');
title('True vs Fused Vertical Position');
legend('True \eta(t)', 'Kalman fused estimate', 'GNSS samples', 'Location', 'best');
xlim([0 120]);

%% Plot fused vertical velocity - first 120 s
figure;
plot(t, v, 'LineWidth', 1.2); hold on;
plot(t, v_kf, 'LineWidth', 1.1);
grid on;
xlabel('Time [s]');
ylabel('Vertical velocity [m/s]');
title('True vs Fused Vertical Velocity');
legend('True v(t)', 'Kalman fused estimate', 'Location', 'best');
xlim([0 120]);

%% Zoomed plots - first 20 s
figure;
subplot(2,1,1);
plot(t(t<=20), eta(t<=20), 'LineWidth', 1.2); hold on;
plot(t(t<=20), eta_kf(t<=20), 'LineWidth', 1.1);
plot(t_gnss(t_gnss<=20), eta_gnss_meas(t_gnss<=20), 'o');
grid on;
ylabel('\eta [m]');
title('Fused Vertical Position - First 20 s');
legend('True', 'Fused', 'GNSS', 'Location', 'best');

subplot(2,1,2);
plot(t(t<=20), v(t<=20), 'LineWidth', 1.2); hold on;
plot(t(t<=20), v_kf(t<=20), 'LineWidth', 1.1);
grid on;
xlabel('Time [s]');
ylabel('v [m/s]');
title('Fused Vertical Velocity - First 20 s');
legend('True', 'Fused', 'Location', 'best');

%% Fusion error statistics
eta_err_kf = eta_kf - eta;
v_err_kf   = v_kf - v;

fprintf('\n--- 1D IMU-GNSS fusion results ---\n');
fprintf('Position mean error: %.6e m\n', mean(eta_err_kf));
fprintf('Position std error:  %.6e m\n', std(eta_err_kf));
fprintf('Position max abs error: %.6e m\n', max(abs(eta_err_kf)));

fprintf('Velocity mean error: %.6e m/s\n', mean(v_err_kf));
fprintf('Velocity std error:  %.6e m/s\n', std(v_err_kf));
fprintf('Velocity max abs error: %.6e m/s\n', max(abs(v_err_kf)));

return