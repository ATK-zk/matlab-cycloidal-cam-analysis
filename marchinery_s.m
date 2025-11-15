%% ======================== CAM Profile Analysis: Cycloidal Motion ========================
% วิเคราะห์โปรไฟล์ cam แบบ cycloidal สำหรับระบบวาล์ว Low-Lift & High-Lift

clear; clc; format long g;
set(0, 'DefaultAxesFontSize', 11, 'DefaultAxesFontName', 'Arial', ...
       'DefaultTextFontSize', 11, 'DefaultTextFontName', 'Arial');

%% ========================= ค่าคงที่และพารามิเตอร์ =============================

deg = pi/180;               % Conversion factor: degree to radian
L_low  = 10e-3;             % Valve lift Low-Lift [m]
L_high = 10e-3;             % Valve lift High-Lift [m]

% Low-Lift profile angles
theta0_low      = 100 * deg;
beta_low_rise   = 52.75 * deg;
beta_low_return = 52.75 * deg;
theta1_low      = theta0_low + beta_low_rise;
theta2_low      = theta1_low;
theta3_low      = theta2_low + beta_low_return;

% High-Lift profile angles
theta0_high     = 130 * deg;
beta_high_rise  = 64.6 * deg;
beta_high_return = 64.6 * deg;
theta1_high     = theta0_high + beta_high_rise;
theta2_high     = theta1_high;
theta3_high     = theta2_high + beta_high_return;

%% ======================== Symbolic Definitions & Derivatives ==============

syms theta omega real

% Cycloidal motion functions
rise = @(th, th0, beta, L) L * ((th - th0) / beta - sin(2*pi*(th - th0)/beta) / (2*pi));
ret  = @(th, th2, beta, L) L * (1 - (th - th2) / beta + sin(2*pi*(th - th2)/beta) / (2*pi));

% Position profiles (symbolic piecewise)
s_low = piecewise(...
    theta < theta0_low,        0, ...
    theta < theta1_low,        rise(theta, theta0_low, beta_low_rise, L_low), ...
    theta < theta3_low,        ret(theta, theta2_low, beta_low_return, L_low), ...
    theta < 2*pi,              0, 0);

s_high = piecewise(...
    theta < theta0_high,       0, ...
    theta < theta1_high,       rise(theta, theta0_high, beta_high_rise, L_high), ...
    theta < theta3_high,       ret(theta, theta2_high, beta_high_return, L_high), ...
    theta < 2*pi,              0, 0);

% Derivatives
ds_dtheta_low   = diff(s_low, theta);
d2s_dtheta_low  = diff(s_low, theta, 2);
ds_dtheta_high  = diff(s_high, theta);
d2s_dtheta_high = diff(s_high, theta, 2);

% Convert to time derivatives: v = (ds/dθ)*ω, a = (d²s/dθ²)*ω²
v_low   = ds_dtheta_low  * omega;
a_low   = d2s_dtheta_low * omega^2;
v_high  = ds_dtheta_high * omega;
a_high  = d2s_dtheta_high * omega^2;

% Generate numeric functions
matlabFunction(s_low,  'Vars', {theta},        'File', 's_low_func');
matlabFunction(s_high, 'Vars', {theta},        'File', 's_high_func');
matlabFunction(v_low,  'Vars', {theta, omega}, 'File', 'v_low_func');
matlabFunction(v_high, 'Vars', {theta, omega}, 'File', 'v_high_func');
matlabFunction(a_low,  'Vars', {theta, omega}, 'File', 'a_low_func');
matlabFunction(a_high, 'Vars', {theta, omega}, 'File', 'a_high_func');

% Create function handles with vectorization support
wrapTheta = @(th) mod(th, 2*pi);
s_low_fun  = @(th)     arrayfun(@(x) s_low_func( wrapTheta(x) ), th);
v_low_fun  = @(th, om) arrayfun(@(x) v_low_func( wrapTheta(x), om ), th);
a_low_fun  = @(th, om) arrayfun(@(x) a_low_func( wrapTheta(x), om ), th);
s_high_fun = @(th)     arrayfun(@(x) s_high_func( wrapTheta(x) ), th);
v_high_fun = @(th, om) arrayfun(@(x) v_high_func( wrapTheta(x), om ), th);
a_high_fun = @(th, om) arrayfun(@(x) a_high_func( wrapTheta(x), om ), th);

%% ====================== Data Preparation at 1000 rpm ======================

omega_val = 2*pi * 1000 / 60;       % Convert 1000 rpm to rad/s
theta_samp = linspace(0, 2*pi, 1500);
theta_deg = theta_samp * 180 / pi;

% Calculate follower motion
s_low_num  = s_low_fun(theta_samp);
v_low_num  = v_low_fun(theta_samp, omega_val);
a_low_num  = a_low_fun(theta_samp, omega_val);
s_high_num = s_high_fun(theta_samp);
v_high_num = v_high_fun(theta_samp, omega_val);
a_high_num = a_high_fun(theta_samp, omega_val);

% Lever arm ratio (follower to valve)
ratio = 37.4 / 38.3;
s_valve_low  = -ratio * s_low_num;
v_valve_low  = -ratio * v_low_num;
a_valve_low  = -ratio * a_low_num;
s_valve_high = -ratio * s_high_num;
v_valve_high = -ratio * v_high_num;
a_valve_high = -ratio * a_high_num;

%% ======================== Plotting (3x2 subplots) ==========================

fig = figure('Position', [100, 100, 1200, 850]);

% Follower - Position
subplot(3,2,1);
plot(theta_deg, s_low_num*1e3, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_deg, s_high_num*1e3, 'r-', 'LineWidth', 1.5);
ylabel('Lift (mm)', 'FontWeight', 'bold');
title('Follower Displacement', 'FontWeight', 'bold');
legend('Low', 'High', 'Location', 'Best', 'FontSize', 10);
grid on; xlim([0 360]);

% Follower - Velocity
subplot(3,2,2);
plot(theta_deg, v_low_num, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_deg, v_high_num, 'r-', 'LineWidth', 1.5);
ylabel('Velocity (m/s)', 'FontWeight', 'bold');
title('Follower Velocity', 'FontWeight', 'bold');
legend('Low', 'High', 'Location', 'Best', 'FontSize', 10);
grid on; xlim([0 360]);

% Follower - Acceleration
subplot(3,2,3);
plot(theta_deg, a_low_num, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_deg, a_high_num, 'r-', 'LineWidth', 1.5);
ylabel('Acceleration (m/s²)', 'FontWeight', 'bold');
title('Follower Acceleration', 'FontWeight', 'bold');
legend('Low', 'High', 'Location', 'Best', 'FontSize', 10);
grid on; xlim([0 360]);
xlabel('Cam Angle (°)', 'FontWeight', 'bold');

% Valve - Position
subplot(3,2,4);
plot(theta_deg, s_valve_low*1e3, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_deg, s_valve_high*1e3, 'r-', 'LineWidth', 1.5);
ylabel('Lift (mm)', 'FontWeight', 'bold');
title('Valve Displacement', 'FontWeight', 'bold');
legend('Low', 'High', 'Location', 'Best', 'FontSize', 10);
grid on; xlim([0 360]);

% Valve - Velocity
subplot(3,2,5);
plot(theta_deg, v_valve_low, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_deg, v_valve_high, 'r-', 'LineWidth', 1.5);
ylabel('Velocity (m/s)', 'FontWeight', 'bold');
title('Valve Velocity', 'FontWeight', 'bold');
legend('Low', 'High', 'Location', 'Best', 'FontSize', 10);
grid on; xlim([0 360]);
xlabel('Cam Angle (°)', 'FontWeight', 'bold');

% Valve - Acceleration
subplot(3,2,6);
plot(theta_deg, a_valve_low, 'b-', 'LineWidth', 1.5); hold on;
plot(theta_deg, a_valve_high, 'r-', 'LineWidth', 1.5);
ylabel('Acceleration (m/s²)', 'FontWeight', 'bold');
title('Valve Acceleration', 'FontWeight', 'bold');
legend('Low', 'High', 'Location', 'Best', 'FontSize', 10);
grid on; xlim([0 360]);
xlabel('Cam Angle (°)', 'FontWeight', 'bold');

%% ======================== Summary Results ==================================

% Find peak values and their positions
[max_s_low,    idx_s_low]    = max(s_low_num);
[max_s_high,   idx_s_high]   = max(s_high_num);
[max_v_low,    idx_v_low]    = max(abs(v_low_num));
[max_v_high,   idx_v_high]   = max(abs(v_high_num));
[max_a_low,    idx_a_low]    = max(abs(a_low_num));
[max_a_high,   idx_a_high]   = max(abs(a_high_num));

[max_sv_low,   idx_sv_low]   = max(abs(s_valve_low));
[max_sv_high,  idx_sv_high]  = max(abs(s_valve_high));
[max_vv_low,   idx_vv_low]   = max(abs(v_valve_low));
[max_vv_high,  idx_vv_high]  = max(abs(v_valve_high));
[max_av_low,   idx_av_low]   = max(abs(a_valve_low));
[max_av_high,  idx_av_high]  = max(abs(a_valve_high));

% Corresponding angles
angle_mat = @(idx) theta_deg(idx);

fprintf('\n%s\n', repmat('=', 1, 65));
fprintf('  CAM PROFILE ANALYSIS @ 1000 RPM\n');
fprintf('%s\n', repmat('=', 1, 65));

fprintf('\n  FOLLOWER MOTION (Low-Lift start @ 100° High-Lift start @ 130°) \n');
fprintf('  %s\n', repmat('-', 1, 61));
fprintf('  %-20s  %12s  %12s  %10s\n', 'Parameter', 'Low-Lift', 'High-Lift', 'Angle (°)');
fprintf('  %s\n', repmat('-', 1, 61));
fprintf('  Max Displacement    %10.3f mm  %10.3f mm   [%.1f° / %.1f°]\n', ...
    max_s_low*1e3, max_s_high*1e3, angle_mat(idx_s_low), angle_mat(idx_s_high));
fprintf('  Max Velocity        %10.4f m/s %10.4f m/s  [%.1f° / %.1f°]\n', ...
    max_v_low, max_v_high, angle_mat(idx_v_low), angle_mat(idx_v_high));
fprintf('  Max Acceleration    %10.1f m/s² %10.1f m/s² [%.1f° / %.1f°]\n', ...
    max_a_low, max_a_high, angle_mat(idx_a_low), angle_mat(idx_a_high));

fprintf('\n  VALVE MOTION \n');
fprintf('  %s\n', repmat('-', 1, 61));
fprintf('  Max Displacement    %10.3f mm  %10.3f mm   [%.1f° / %.1f°]\n', ...
    max_sv_low*1e3, max_sv_high*1e3, angle_mat(idx_sv_low), angle_mat(idx_sv_high));
fprintf('  Max Velocity        %10.4f m/s %10.4f m/s  [%.1f° / %.1f°]\n', ...
    max_vv_low, max_vv_high, angle_mat(idx_vv_low), angle_mat(idx_vv_high));
fprintf('  Max Acceleration    %10.1f m/s² %10.1f m/s² [%.1f° / %.1f°]\n', ...
    max_av_low, max_av_high, angle_mat(idx_av_low), angle_mat(idx_av_high));
fprintf('\n%s\n\n', repmat('=', 1, 65));
