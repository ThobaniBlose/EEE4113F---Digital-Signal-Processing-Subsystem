%% Heave Signal Generation: Displacement, Velocity, and Acceleration
% This script generates a clean synthetic heave displacement signal:
%
%   eta(t) = sum_{i=1}^{N} A_i sin(2*pi*f_i*t + phi_i)
%
% and then computes the corresponding true velocity and acceleration:
%
%   v(t) = d(eta)/dt
%   a(t) = d^2(eta)/dt^2
%
% The goal here is not realism yet.
% The goal is to create a known ground-truth signal for later PSD
% and wave-parameter extraction.

clc;
clear;
close all;

%% 1. Simulation settings
% Create a time vector for 30 minutes sampled at 100 Hz.

fs = 100;                      % Sampling frequency [Hz]
dt = 1/fs;                     % Sampling interval [s]
duration_min = 30;             % Duration [minutes]
duration_s = duration_min*60;  % Duration [s]

t = 0:dt:(duration_s - dt);    % Time vector (exactly 30 minutes)

%% 2. Define wave components
% Choose a few simple sinusoidal components in a rough ocean-wave range.
% These values are kept small and simple on purpose.

N = 4;   % Number of sinusoidal components

A   = [0.08 0.05 0.035 0.02];      % Amplitudes [m]
f   = [0.07 0.11 0.16 0.22];       % Frequencies [Hz]
phi = [0 pi/6 pi/3 pi/2];          % Phases [rad]

%% 3. Initialise signals
% Preallocate arrays for displacement, velocity, and acceleration.

eta = zeros(size(t));   % Displacement [m]
v   = zeros(size(t));   % Velocity [m/s]
a   = zeros(size(t));   % Acceleration [m/s^2]

%% 4. Generate displacement, true velocity, and true acceleration
% For each sinusoidal component:
%
%   eta_i(t) = A_i sin(2*pi*f_i*t + phi_i)
%
% Differentiate analytically:
%
%   v_i(t) = A_i (2*pi*f_i) cos(2*pi*f_i*t + phi_i)
%
%   a_i(t) = -A_i (2*pi*f_i)^2 sin(2*pi*f_i*t + phi_i)

for i = 1:N
    omega_i = 2*pi*f(i);   % Angular frequency [rad/s]
    
    % Displacement contribution
    eta = eta + A(i)*sin(omega_i*t + phi(i));
    
    % Velocity contribution
    v = v + A(i)*omega_i*cos(omega_i*t + phi(i));
    
    % Acceleration contribution
    a = a - A(i)*(omega_i^2)*sin(omega_i*t + phi(i));
end

%% 5. Plot displacement, velocity, and acceleration
% Plot all three true signals.

figure;

subplot(3,1,1);
plot(t, eta, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('\eta(t) [m]');
title('True Heave Displacement');

subplot(3,1,2);
plot(t, v, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('v(t) [m/s]');
title('True Heave Velocity');

subplot(3,1,3);
plot(t, a, 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('a(t) [m/s^2]');
title('True Heave Acceleration');

%% 6. Optional: zoom into a short section for easier viewing
% Over 30 minutes the full signal is long, so it is useful to inspect
% just a short segment as well.

zoom_duration = 120;                         % Duration to view [s]
idx_zoom = t <= zoom_duration;              % Indices for first 120 s

figure;

subplot(3,1,1);
plot(t(idx_zoom), eta(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('\eta(t) [m]');
title('True Heave Displacement (First 120 s)');

subplot(3,1,2);
plot(t(idx_zoom), v(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('v(t) [m/s]');
title('True Heave Velocity (First 120 s)');

subplot(3,1,3);
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('a(t) [m/s^2]');
title('True Heave Acceleration (First 120 s)');

%% 7. Display chosen component parameters
% Show the amplitudes, frequencies, and phases used.

disp('Chosen sinusoidal wave components:');
disp(table((1:N)', A', f', phi', ...
    'VariableNames', {'Component', 'Amplitude_m', 'Frequency_Hz', 'Phase_rad'}));

%% 8. Ground-truth PSD from true displacement eta(t)
% Use Welch PSD on the known displacement signal

window_length = 4096;                 % choose a reasonable segment length
overlap = window_length/2;
nfft = 4096;

[PSD_eta, f_psd] = pwelch(eta, hann(window_length), overlap, nfft, fs);

%% 9. Restrict to wave band
% Use a wave band consistent with your synthetic components

f_low = 0.03;
f_high = 0.40;

idx_band = (f_psd >= f_low) & (f_psd <= f_high);

f_wave = f_psd(idx_band);
S_wave = PSD_eta(idx_band);

%% 10. Spectral moments
m0 = trapz(f_wave, S_wave);
m1 = trapz(f_wave, f_wave .* S_wave);
m2 = trapz(f_wave, (f_wave.^2) .* S_wave);

%% 11. Bulk wave parameters
Hm0  = 4*sqrt(m0);
Tm01 = m0/m1;
T0   = sqrt(m0/m2);

[~, idx_peak] = max(S_wave);
fp = f_wave(idx_peak);
Tp = 1/fp;

%% 12. Plot PSD
figure;
plot(f_psd, PSD_eta, 'LineWidth', 1.2);
grid on;
xlabel('Frequency [Hz]');
ylabel('S_{\eta}(f) [m^2/Hz]');
title('Ground-Truth Displacement PSD from \eta(t)');

xlim([0 0.5]);

%% 13. Display results
fprintf('\nGround-truth wave parameters from eta(t):\n');
fprintf('Hm0  = %.4f m\n', Hm0);
fprintf('Tm01 = %.4f s\n', Tm01);
fprintf('T0   = %.4f s\n', T0);
fprintf('Tp   = %.4f s\n', Tp);
fprintf('fp   = %.4f Hz\n', fp);

%% 14. Simulated IMU acceleration measurement
% Create a simple measured acceleration:
% a_meas(t) = a_true(t) + bias + noise

bias = 0.005;              % constant bias [m/s^2]
noise_std = 0.01;          % white noise std [m/s^2]

rng(1);                    % for repeatability
noise = noise_std * randn(size(a));

a_meas = a + bias + noise;

%% 15. Plot comparison
zoom_duration = 60;                    % first 60 s
idx_zoom = t <= zoom_duration;

figure;

subplot(2,1,1);
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('a(t) [m/s^2]');
title('True Heave Acceleration');

subplot(2,1,2);
plot(t(idx_zoom), a_meas(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('a_{meas}(t) [m/s^2]');
title('Simulated IMU Measured Acceleration');

%% 16. Show bias/noise values
fprintf('\nSimulated IMU measurement model:\n');
fprintf('Bias      = %.4f m/s^2\n', bias);
fprintf('Noise std = %.4f m/s^2\n', noise_std);

%% 17. Bias removal from measured acceleration
% Remove the DC component / constant bias estimate

bias_est = mean(a_meas);
a_corr = a_meas - bias_est;

%% 18. Plot corrected acceleration
zoom_duration = 60;
idx_zoom = t <= zoom_duration;

figure;

subplot(3,1,1);
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('a(t) [m/s^2]');
title('True Heave Acceleration');

subplot(3,1,2);
plot(t(idx_zoom), a_meas(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('a_{meas}(t) [m/s^2]');
title('Measured Acceleration');

subplot(3,1,3);
plot(t(idx_zoom), a_corr(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('a_{corr}(t) [m/s^2]');
title('Bias-Corrected Acceleration');

%% 19. Print bias comparison
fprintf('\nBias handling:\n');
fprintf('True bias added      = %.6f m/s^2\n', bias);
fprintf('Estimated bias mean  = %.6f m/s^2\n', bias_est);
fprintf('Bias estimation error = %.6f m/s^2\n', bias_est - bias);

%% 20. Band-pass filter the bias-corrected acceleration
% Keep the wave band and suppress very low-frequency drift
% and high-frequency noise.
%
% Use designfilt instead of butter(...,'sos') for compatibility.

f_low  = 0.03;    % lower cutoff [Hz]
f_high = 0.40;    % upper cutoff [Hz]

d = designfilt('bandpassiir', ...
    'FilterOrder', 4, ...
    'HalfPowerFrequency1', f_low, ...
    'HalfPowerFrequency2', f_high, ...
    'SampleRate', fs);

% Zero-phase filtering
a_filt = filtfilt(d, a_corr);

%% 21. Check filtered signal is valid
fprintf('\nFilter validity checks:\n');
fprintf('Any NaNs in a_filt? %d\n', any(isnan(a_filt)));
fprintf('Any Infs in a_filt? %d\n', any(isinf(a_filt)));

%% 22. Plot filtered acceleration
zoom_duration = 60;
idx_zoom = t <= zoom_duration;

figure;

subplot(3,1,1);
plot(t(idx_zoom), a(idx_zoom), 'LineWidth', 1.2);
grid on;
xlabel('Time [s]');
ylabel('a(t) [m/s^2]');
title('True Heave Acceleration');

subplot(3,1,2);
plot(t(idx_zoom), a_corr(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('a_{corr}(t) [m/s^2]');
title('Bias-Corrected Acceleration');

subplot(3,1,3);
plot(t(idx_zoom), a_filt(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('a_{filt}(t) [m/s^2]');
title('Band-Pass Filtered Acceleration');

%% 23. PSD comparison before and after filtering
[PSD_corr, f_corr] = pwelch(a_corr, hann(4096), 2048, 4096, fs);
[PSD_filt, f_filt] = pwelch(a_filt, hann(4096), 2048, 4096, fs);

figure;
plot(f_corr, PSD_corr, 'LineWidth', 1.0);
hold on;
plot(f_filt, PSD_filt, 'LineWidth', 1.0);
grid on;
xlabel('Frequency [Hz]');
ylabel('S_a(f) [(m/s^2)^2/Hz]');
title('PSD of Acceleration Before and After Filtering');
legend('Bias-corrected', 'Filtered');
xlim([0 1]);

%% 24. First integration: acceleration to velocity
v_est = cumtrapz(t, a_filt);

%% 25. Detrend estimated velocity
% Remove constant + linear trend after integration
v_est = detrend(v_est, 'linear');

%% 26. Compare away from startup transient
t_start_compare = 200;              % ignore first 200 s
idx_cmp = t >= t_start_compare & t <= (t_start_compare + 60);

figure;
plot(t(idx_cmp), v(idx_cmp), 'LineWidth', 1.2); hold on;
plot(t(idx_cmp), v_est(idx_cmp), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('True vs Estimated Velocity (after startup transient)');
legend('True velocity', 'Estimated velocity');
%% 27. Compare true and estimated velocity on same plot
figure;
plot(t(idx_zoom), v(idx_zoom), 'LineWidth', 1.2); hold on;
plot(t(idx_zoom), v_est(idx_zoom), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Velocity [m/s]');
title('True vs Estimated Velocity');
legend('True velocity', 'Estimated velocity');

%% 28. Second integration: velocity to displacement
eta_est = cumtrapz(t, v_est);

%% 29. Detrend estimated displacement
eta_est = detrend(eta_est, 'linear');

%% 29.5 Remove residual displacement drift
% After two integrations, residual low-frequency drift can remain.
% Remove it with linear detrending.

eta_est_corr = detrend(eta_est, 'linear');

%% 30. Compare true vs corrected estimated displacement
t_start_compare = 200;              % ignore first 200 s
idx_cmp = t >= t_start_compare & t <= (t_start_compare + 60);

figure;
plot(t(idx_cmp), eta(idx_cmp), 'LineWidth', 1.2); hold on;
plot(t(idx_cmp), eta_est_corr(idx_cmp), 'LineWidth', 1.0);
grid on;
xlabel('Time [s]');
ylabel('Displacement [m]');
title('True vs Estimated Displacement (detrended)');
legend('True displacement', 'Estimated displacement');

%% 31. Prepare estimated displacement for spectral analysis
% Remove mean before PSD estimation

eta_est_psd = eta_est_corr - mean(eta_est_corr);

%% 32. PSD of estimated displacement
window_length = 4096;
overlap = window_length/2;
nfft = 4096;

[PSD_eta_est, f_eta_est] = pwelch(eta_est_psd, hann(window_length), overlap, nfft, fs);

%% 33. Plot ground-truth vs estimated displacement PSD
figure;
plot(f_psd, PSD_eta, 'LineWidth', 1.2); hold on;
plot(f_eta_est, PSD_eta_est, 'LineWidth', 1.0);
grid on;
xlabel('Frequency [Hz]');
ylabel('S_{\eta}(f) [m^2/Hz]');
title('Ground-Truth vs Estimated Displacement PSD');
legend('Ground-truth PSD', 'Estimated PSD');
xlim([0 0.5]);

%% 34. Restrict estimated PSD to wave band
idx_band_est = (f_eta_est >= f_low) & (f_eta_est <= f_high);

f_wave_est = f_eta_est(idx_band_est);
S_wave_est = PSD_eta_est(idx_band_est);

%% 35. Spectral moments from estimated displacement
m0_est = trapz(f_wave_est, S_wave_est);
m1_est = trapz(f_wave_est, f_wave_est .* S_wave_est);
m2_est = trapz(f_wave_est, (f_wave_est.^2) .* S_wave_est);

%% 36. Bulk wave parameters from estimated displacement
Hm0_est  = 4*sqrt(m0_est);
Tm01_est = m0_est/m1_est;
T0_est   = sqrt(m0_est/m2_est);

[~, idx_peak_est] = max(S_wave_est);
fp_est = f_wave_est(idx_peak_est);
Tp_est = 1/fp_est;

%% 37. Display comparison
fprintf('\nEstimated wave parameters from reconstructed displacement:\n');
fprintf('Hm0_est  = %.4f m\n', Hm0_est);
fprintf('Tm01_est = %.4f s\n', Tm01_est);
fprintf('T0_est   = %.4f s\n', T0_est);
fprintf('Tp_est   = %.4f s\n', Tp_est);
fprintf('fp_est   = %.4f Hz\n', fp_est);

fprintf('\nGround-truth vs estimated:\n');
fprintf('Hm0  : %.4f -> %.4f m\n', Hm0, Hm0_est);
fprintf('Tm01 : %.4f -> %.4f s\n', Tm01, Tm01_est);
fprintf('T0   : %.4f -> %.4f s\n', T0, T0_est);
fprintf('Tp   : %.4f -> %.4f s\n', Tp, Tp_est);

%% 38. Better PSD comparison in the wave band
figure;
plot(f_psd(idx_band), PSD_eta(idx_band), 'LineWidth', 1.2); hold on;
plot(f_wave_est, S_wave_est, 'LineWidth', 1.0);
grid on;
xlabel('Frequency [Hz]');
ylabel('S_{\eta}(f) [m^2/Hz]');
title('Ground-Truth vs Estimated Displacement PSD (Wave Band Only)');
legend('Ground-truth PSD', 'Estimated PSD');
xlim([0.03 0.30]);