%% Two Independent Cam Lobes (Separate Profiles, Shared Base Circle)
clear; clc; close all; format long g;

%% -------------------- Shared base circle --------------------
R_base    = 20e-3;                     % Base circle radius
theta_samp = linspace(0, 2*pi, 1200);
x_base = R_base * cos(theta_samp);
y_base = R_base * sin(theta_samp);

L_tan = 30e-3;                         % Tangent length for construction
phi    = linspace(0, 2*pi, 2400);      % Ellipse sampling

%% -------------------- Lobe 1 parameters --------------------
theta_start1 = (180+64.6)*pi/180;
theta_end1   = (180-64.6)*pi/180;
a1 = 10e-3; b1 = 5e-3; offset1 = 25e-3; theta_peak1 = 180*pi/180;

%% -------------------- Lobe 2 parameters --------------------
theta_start2 = (180+53.1)*pi/180;
theta_end2   = (180-53.1)*pi/180;
a2 = 5e-3; b2 = 5e-3; offset2 = 25e-3; theta_peak2 = 180*pi/180;
%% -------------------- Build Lobe 1 --------------------
[x_lobe1, y_lobe1] = build_lobe(theta_start1, theta_end1, R_base, L_tan, ...
    x_base, y_base, a1, b1, offset1, theta_peak1, phi);

%% -------------------- Build Lobe 2 --------------------
[x_lobe2, y_lobe2] = build_lobe(theta_start2, theta_end2, R_base, L_tan, ...
    x_base, y_base, a2, b2, offset2, theta_peak2, phi);

%% -------------------- Plot lobes separately --------------------
figure('Name','Independent Cam Lobes'); hold on;
plot(x_lobe1*1e3, y_lobe1*1e3, 'g-', 'LineWidth', 2);
plot(x_lobe2*1e3, y_lobe2*1e3, 'c-', 'LineWidth', 2);
axis equal; grid on;
xlabel('X (mm)'); ylabel('Y (mm)');
title('Two Independent Cam Lobes (Separate Profiles)');
legend('Lobe 1','Lobe 2');

%% ==================== Local function ====================
function [x_lobe, y_lobe] = build_lobe(theta_start, theta_end, R_base, L_tan, ...
    x_base, y_base, a_ellipse, b_ellipse, offset, theta_peak, phi)

    % Ellipse construction
    ellipse_center = offset * [cos(theta_peak); sin(theta_peak)];
    theta_rot = theta_peak - pi/2;
    Rrot = [cos(theta_rot), -sin(theta_rot); sin(theta_rot), cos(theta_rot)];
    x_local = a_ellipse * cos(phi);
    y_local = b_ellipse * sin(phi);
    E = Rrot * [x_local; y_local];
    x_ellipse = E(1,:) + ellipse_center(1);
    y_ellipse = E(2,:) + ellipse_center(2);

    % Tangent points on base circle
    x0_start = R_base * cos(theta_start); y0_start = R_base * sin(theta_start);
    x0_end   = R_base * cos(theta_end);   y0_end   = R_base * sin(theta_end);

    % Tangent directions
    dx_start = -sin(theta_start); dy_start = cos(theta_start);
    dx_end   = -sin(theta_end);   dy_end   = cos(theta_end);
    t1_dir = [dx_start; dy_start] / hypot(dx_start, dy_start);
    t2_dir = [dx_end;   dy_end]   / hypot(dx_end,   dy_end);

    % Intersections with ellipse
    tan1_pt = [x0_start; y0_start];
    tan2_pt = [x0_end;   y0_end];
    rel1 = [x_ellipse; y_ellipse] - tan1_pt;
    rel2 = [x_ellipse; y_ellipse] - tan2_pt;
    dist1 = rel1(1,:).*t1_dir(2) - rel1(2,:).*t1_dir(1);
    dist2 = rel2(1,:).*t2_dir(2) - rel2(2,:).*t2_dir(1);

    c1 = find(abs(diff(sign(dist1))) > 1);
    c2 = find(abs(diff(sign(dist2))) > 1);
    if isempty(c1) || isempty(c2)
        error('Ellipse does not intersect both tangent lines for given angles.');
    end
    i1 = c1(1); i2 = c2(1);

    % Refine intersection points
    tlin1 = dist1(i1) / (dist1(i1) - dist1(i1+1));
    p1 = [x_ellipse(i1) + tlin1*(x_ellipse(i1+1)-x_ellipse(i1));
          y_ellipse(i1) + tlin1*(y_ellipse(i1+1)-y_ellipse(i1))];
    tlin2 = dist2(i2) / (dist2(i2) - dist2(i2+1));
    p2 = [x_ellipse(i2) + tlin2*(x_ellipse(i2+1)-x_ellipse(i2));
          y_ellipse(i2) + tlin2*(y_ellipse(i2+1)-y_ellipse(i2))];

    % Arc selection: take the ellipse arc that lies outside the base circle
    r_ellipse = hypot(x_ellipse, y_ellipse);
    outside_mask = r_ellipse > R_base;
    idx_candidates = find(outside_mask);
    % keep only between intersections
    if i1 < i2
        idx_arc = i1:i2;
    else
        idx_arc = i2:i1;
    end
    idx_arc = intersect(idx_arc, idx_candidates);
    x_arc = x_ellipse(idx_arc);
    y_arc = y_ellipse(idx_arc);

    % Tangent segments
    x_t1 = linspace(x0_start, p1(1), 80);
    y_t1 = linspace(y0_start,  p1(2), 80);
    x_t2 = linspace(p2(1), x0_end,  80);
    y_t2 = linspace(p2(2), y0_end,   80);

    % Trim base circle
    ang_base = mod(atan2(y_base, x_base), 2*pi);
    theta1 = mod(theta_start, 2*pi);
    theta2 = mod(theta_end,   2*pi);
    delta = mod(theta2 - theta1, 2*pi);
    if delta < pi
        between = mod(ang_base - theta1, 2*pi) <= delta;
    else
        span = 2*pi - delta;
        between = mod(ang_base - theta2, 2*pi) <= span;
    end
    mask_keep_base = ~between;
    x_base_trim = x_base(mask_keep_base);
    y_base_trim = y_base(mask_keep_base);

    % Assemble lobe: base → tan1 → arc → tan2
    x_lobe = [x_base_trim, x_t2(2:end), x_arc(2:end), x_t1(2:end)];
    y_lobe = [y_base_trim, y_t2(2:end), y_arc(2:end), y_t1(2:end)];
end


%% Lift vs Cam Angle (0–360°) for Two Independent Lobes
% Assumes x_lobe1, y_lobe1, x_lobe2, y_lobe2 exist from prior construction.

R_base = 20e-3;  % base circle radius (same as used in lobe scripts)

% --- Helper: compute angle (0–360°) and lift (mm), sorted by angle ---
function [theta_deg_sorted, lift_mm_sorted] = lift_profile_0_360(x_lobe, y_lobe, R_base)
    theta_deg = mod(rad2deg(atan2(y_lobe, x_lobe)), 360);  % wrap to [0,360)
    r = hypot(x_lobe, y_lobe);
    lift_mm = (r - R_base) * 1e3;                          % mm
    [theta_deg_sorted, idx] = sort(theta_deg);
    lift_mm_sorted = lift_mm(idx);
end

% --- Lobe 1 ---
[theta1_deg, lift1_mm] = lift_profile_0_360(x_lobe1, y_lobe1, R_base);

% --- Lobe 2 ---
[theta2_deg, lift2_mm] = lift_profile_0_360(x_lobe2, y_lobe2, R_base);

% --- Plot both (0–360°) ---
figure('Name','Lift vs Cam Angle (0–360°)'); hold on;
plot(theta1_deg, lift1_mm, 'b-', 'LineWidth', 2);
plot(theta2_deg, lift2_mm, 'r-', 'LineWidth', 2);
grid on; xlim([0 360]);
xlabel('Cam Angle (deg)'); ylabel('Lift (mm)');
title('Lift vs Cam Angle (0–360°) for Two Independent Lobes');
legend('Lobe 1','Lobe 2','Location','best');

%% ==================== Lift Analysis and Dynamics ====================
% Continue from your lobe construction + lift plotting code

omega_val = 2*pi*1000/60;   % rad/s at 1000 rpm
ratio = 37.4/38.3;          % lever arm ratio follower->valve

% Interpolate lift onto uniform angle grid [0,360]
theta_grid = linspace(0,360,1500);
theta_rad  = theta_grid*pi/180;

% Interpolation helper
lift_interp = @(theta_deg,lift_mm,theta_grid) ...
    interp1(theta_deg, lift_mm/1e3, theta_grid, 'linear', 0); % convert mm->m

lift1 = lift_interp(theta1_deg,lift1_mm,theta_grid); % lobe1 lift [m]
lift2 = lift_interp(theta2_deg,lift2_mm,theta_grid); % lobe2 lift [m]

% Time step per sample
dt = (theta_rad(2)-theta_rad(1))/omega_val;

% Follower velocity/acceleration
v1 = gradient(lift1,dt); a1 = gradient(v1,dt);
v2 = gradient(lift2,dt); a2 = gradient(v2,dt);

% Valve motion
s_valve1 = -ratio*lift1; v_valve1 = -ratio*v1; a_valve1 = -ratio*a1;
s_valve2 = -ratio*lift2; v_valve2 = -ratio*v2; a_valve2 = -ratio*a2;

%% ==================== Plotting ====================

% Follower motion (Position -> Velocity -> Acceleration)
figure('Position',[100,100,800,900]);

subplot(3,1,1);
plot(theta_grid,lift2*1e3,'b-',theta_grid,lift1*1e3,'r-','LineWidth',1.5);
ylabel('Lift (mm)'); title('Follower Displacement');
legend('Lobe 2','Lobe 1'); grid on; xlim([0 360]);

subplot(3,1,2);
plot(theta_grid,v2,'b-',theta_grid,v1,'r-','LineWidth',1.5);
ylabel('Velocity (m/s)'); title('Follower Velocity');
legend('Lobe 2','Lobe 1'); grid on; xlim([0 360]);

subplot(3,1,3);
plot(theta_grid,a2,'b-',theta_grid,a1,'r-','LineWidth',1.5);
ylabel('Acceleration (m/s²)'); title('Follower Acceleration');
legend('Lobe 2','Lobe 1'); grid on; xlim([0 360]); xlabel('Cam Angle (°)');

% Valve motion (Position -> Velocity -> Acceleration)
figure('Position',[100,100,800,900]);

subplot(3,1,1);
plot(theta_grid,s_valve2*1e3,'b-',theta_grid,s_valve1*1e3,'r-','LineWidth',1.5);
ylabel('Lift (mm)'); title('Valve Displacement');
legend('Lobe 2','Lobe 1'); grid on; xlim([0 360]);

subplot(3,1,2);
plot(theta_grid,v_valve2,'b-',theta_grid,v_valve1,'r-','LineWidth',1.5);
ylabel('Velocity (m/s)'); title('Valve Velocity');
legend('Lobe 2','Lobe 1'); grid on; xlim([0 360]);

subplot(3,1,3);
plot(theta_grid,a_valve2,'b-',theta_grid,a_valve1,'r-','LineWidth',1.5);
ylabel('Acceleration (m/s²)'); title('Valve Acceleration');
legend('Lobe 2','Lobe 1'); grid on; xlim([0 360]); xlabel('Cam Angle (°)');

%% ==================== Summary ====================
[max_lift2,idx_s2] = max(lift2); [max_lift1,idx_s1] = max(lift1);
[max_v2,idx_v2] = max(abs(v2));  [max_v1,idx_v1] = max(abs(v1));
[max_a2,idx_a2] = max(abs(a2));  [max_a1,idx_a1] = max(abs(a1));

fprintf('\nCAM PROFILE ANALYSIS @ 1000 RPM\n');
fprintf('-------------------------------------------\n');
fprintf('Follower Max Lift: Lobe2 %.3f mm @ %.1f°, Lobe1 %.3f mm @ %.1f°\n',...
    max_lift2*1e3,theta_grid(idx_s2),max_lift1*1e3,theta_grid(idx_s1));
fprintf('Follower Max Velocity: Lobe2 %.4f m/s, Lobe1 %.4f m/s\n',max_v2,max_v1);
fprintf('Follower Max Accel:    Lobe2 %.1f m/s², Lobe1 %.1f m/s²\n',max_a2,max_a1);
fprintf('Valve Max Lift: Lobe2 %.3f mm, Lobe1 %.3f mm\n',...
    max(abs(s_valve2))*1e3,max(abs(s_valve1))*1e3);
fprintf('Valve Max Velocity: Lobe2 %.4f m/s, Lobe1 %.4f m/s\n',...
    max(abs(v_valve2)),max(abs(v_valve1)));
fprintf('Valve Max Accel:    Lobe2 %.1f m/s², Lobe1 %.1f m/s²\n',...
    max(abs(a_valve2)),max(abs(a_valve1)));
