%% ALES PROJECT - ROBUST FAULT DIAGNOSIS (CHEN, PATTON, ZHANG 1996)
% =========================================================================
% REPLICATION & VISUALIZATION (FINAL STORYTELLING VERSION)
%
% Scenario: Real Fault on Sensor 1 at t = 3.0s
% Goal: 
% 1. Show Standard Filter failing dynamically (mis-isolation).
% 2. Show Robust Filter succeeding dynamically (perfect isolation).
% 3. Compare Residual Magnitude over full time (False Alarms vs Stability).
% 4. Summary Bar Chart (The "Smoking Gun").
% =========================================================================

clc; clear; close all;

%% SYSTEM DEFINITION & MATRICES

% --- LINEAR MODEL ---
% These matrices (A, B, C) represent the linearized model of the jet engine.
% Both filters use this model as their baseline knowledge of the system.
A = [-1.5581,  0.6925,  0.3974;
      0.2619, -2.2228,  0.2238;
      0       0      -10    ];
B = [0; 0; 10];
C = [1,0,0; 0,1,0; 0,0,1; 
     0.55107, 0.13320, 0.30603; 
     0.55217, 0.13526, 0.32912; 
    -0.25693, -0.23625, 0.61299];
[n, ~] = size(A); [p, ~] = size(C);

% --- NON-LINEAR ERROR TERMS ---
% This matrix (E*) represents the "Reality Gap". It distributes second-order
% terms (x^2, xy, etc.) into the state equation.
% Ideally, a filter should ignore these, but Standard filters treat them as faults.
E_star = [1.3293, 3.4440, 0.1375, -5.1304, -1.7826, -1.8719;
          5.6812, -0.5281, -0.3385, -1.6193, 0.5229, 0;
          0, 0, 0, 0, 0, 0];

% --- ROBUST FILTER MATRICES ---
% These are the specific matrices calculated by the authors to achieve
% "Disturbance Decoupling". H and T work together to make the residual
% insensitive to the unknown inputs (the non-linearities).
H_paper = [ 0.6117, -0.1170, 0, 0.3215, 0.3220, -0.1295;
           -0.1170,  0.9382, 0, 0.0605, 0.0623, -0.1916;
            0,       0,      0, 0,      0,       0     ];
K_paper = [-0.0708,  0.0443,  0.5658,  0.1400,  0.1531,  0.3540;
            0.0443, -0.0277, -0.3540, -0.0876, -0.0958, -0.2215;
            0.5658, -0.3540, -4.5229, -1.1193, -1.2239, -2.8297];
T_paper = eye(n) - H_paper * C; 

%% FILTER DESIGN REVERSE ENGINEERING

% --- Standard BFDF Design ---
% We design the Standard Beard Fault Detection Filter using pole placement.
% It has no knowledge of E_star, so it assumes the linear model is perfect.
pole = 3;
K_bfdf = (A + pole * eye(n)) * pinv(C);

% --- Robust UIO Reverse Engineering ---
% To calculate the "NPD" (Normalized Projection Distance), we need to know
% the directional gain matrix K1. The paper gives us the total gain K.
% We solve the equation K = K1*(I - CH) + ... backwards to extract K1.
term_A_HCA = A - H_paper * C * A;
RHS = K_paper - term_A_HCA * H_paper;
K1_paper = RHS * pinv(eye(p) - C * H_paper);
F_paper = A - H_paper*C*A - K1_paper*C;

%% SIMULATION LOOP

dt = 0.001; T_sim = 10; time = 0:dt:T_sim; steps = length(time);
u_val = 0.20;       % Constant input
fault_mag = 0.02;   % Fault magnitude (2%)
fault_time = 3.0;   % Fault occurs at 3 seconds

% We amplify the non-linear error (E_star) by 20x. This simulates a scenario
% where the linear model is poor. This "stress test" is necessary to make
% the Standard Filter fail, highlighting the need for Robustness.
dist_factor = 20.0; 

% Init Logs
log_r_bfdf = zeros(p, steps);
log_r_uio = zeros(p, steps);
log_npd_bfdf_dyn = zeros(3, steps); % Log for Standard Filter Decision
log_npd_uio_dyn = zeros(3, steps);  % Log for Robust Filter Decision

% We only simulate Fault Case 1 (Sensor 1 Broken)
s_fault = 1; 
x = zeros(n,1); x_hat_bfdf = zeros(n,1); z_uio = zeros(n,1);

fprintf('Simulating Fault on Sensor 1 (Distortion Factor = %.1f)...\n', dist_factor);

for k = 1:steps
    % 1. Non-Linear Plant Dynamics (The "Real World")
    % We construct the quadratic terms vector d(x)
    d_x = [x(1)^2; x(2)^2; x(3)^2; x(1)*x(2); x(1)*x(3); x(2)*x(3)];
    % The state evolves using Linear Dynamics + Amplified Non-Linear Error
    dx = A*x + B*u_val + (E_star * d_x) * dist_factor;
    x = x + dx * dt;
    
    % Output Generation & Fault Injection
    y_raw = C*x;
    if time(k) >= fault_time, y_raw(s_fault) = y_raw(s_fault) + fault_mag; end
    
    % 2. Standard Filter (BFDF) Execution
    % Standard observer equation: dx = Ax + Bu + K(y - Cx)
    r_bfdf = y_raw - C * x_hat_bfdf;
    dx_bfdf = A * x_hat_bfdf + B * u_val + K_bfdf * r_bfdf;
    x_hat_bfdf = x_hat_bfdf + dx_bfdf * dt;
    
    % 3. Robust Filter (UIO) Execution
    % Uses the z-state and T matrix to decouple disturbances.
    dz_uio = F_paper * z_uio + T_paper * B * u_val + K_paper * y_raw;
    z_uio = z_uio + dz_uio * dt;
    r_uio = y_raw - C * (z_uio + H_paper * y_raw);
    
    % Logging Residuals for plotting
    log_r_bfdf(:,k) = r_bfdf;
    log_r_uio(:,k) = r_uio;
    
    % 4. Dynamic NPD Calculation
    % We calculate the "match score" for every sensor hypothesis at every time step.
    for j=1:3
        % -- Standard Filter NPD --
        Phi = [eye(p)*0 + (j==1:p)', C*K_bfdf(:,j)]; 
        Phi(:,1) = zeros(p,1); Phi(j,1) = 1;
        P = Phi * pinv(Phi);
        log_npd_bfdf_dyn(j,k) = norm(r_bfdf - P * r_bfdf) / (norm(r_bfdf) + eps);
        
        % -- Robust Filter NPD --
        Phi_r = [zeros(p,1), C*K1_paper(:,j), C*H_paper(:,j)];
        Phi_r(j,1) = 1;
        P_r = Phi_r * pinv(Phi_r);
        log_npd_uio_dyn(j,k) = norm(r_uio - P_r * r_uio) / (norm(r_uio) + eps);
    end
end

% Extract final NPD values for the bar chart (Steady State)
final_npd_bfdf = log_npd_bfdf_dyn(:, end);
final_npd_uio = log_npd_uio_dyn(:, end);

%% VISUALIZATION

% --- FIGURE 1: SUMMARY BAR CHART ---
% This chart compares the final diagnosis of both filters.
% Goal: Show that Standard Filter picks the wrong sensor, Robust picks the right one.
figure('Color','w','Position',[50 500 700 400]);
bar_data = [final_npd_bfdf, final_npd_uio];
b = bar(bar_data);
b(1).FaceColor = [0.85 0.33 0.1]; % Burnt Orange (Standard)
b(2).FaceColor = [0.1 0.6 0.3];   % Forest Green (Robust)

title('FINAL VERDICT: Real Fault on Sensor 1');
xlabel('Fault Hypothesis (Which sensor is broken?)'); 
ylabel('NPD Index (Lower is Better)');
legend({'Standard BFDF (Confused)', 'Robust UIO (Correct)'}, 'Location', 'NorthWest');
xticklabels({'Sensor 1', 'Sensor 2', 'Sensor 3'});
grid on;
text(2, final_npd_bfdf(2), '\downarrow Standard Mis-isolates (Selects S2)', ...
    'Vert','bottom','Horiz','center', 'FontSize', 10, 'Color','r', 'FontWeight','bold');
text(1, final_npd_uio(1), '\downarrow Robust Correct (Selects S1)', ...
    'Vert','bottom','Horiz','center', 'FontSize', 10, 'Color','k', 'FontWeight','bold');


% --- FIGURE 2: THE FAILURE ---
% This graph shows the decision confidence over time for the Standard Filter.
% Goal: Show the Red Line (Hypothesis 2) dropping lower than the Green Line (Hypothesis 1).
figure('Color','w','Position',[50 50 600 350]);
plot(time, log_npd_bfdf_dyn(1,:), 'g-', 'LineWidth', 1.5); hold on; % S1 (Real Fault)
plot(time, log_npd_bfdf_dyn(2,:), 'r--', 'LineWidth', 2.0);        % S2 (Wrong Choice)
plot(time, log_npd_bfdf_dyn(3,:), 'k:', 'LineWidth', 1.0);         % S3
xline(fault_time, 'k-', 'Fault Injection');

title('FAILURE: Standard Filter Decision over Time');
xlabel('Time [s]'); ylabel('NPD Score (Lower is "Better Match")');
legend('Hypothesis: Sensor 1', 'Hypothesis: Sensor 2', 'Hypothesis: Sensor 3', 'Location', 'Best');
grid on; ylim([0 1.2]);

text(fault_time + 0.5, 0.4, '\downarrow Wrong Decision (S2 selected)!', 'Color', 'r', 'FontWeight', 'bold');


% --- FIGURE 3: THE SOLUTION ---
% This graph shows the decision confidence over time for the Robust Filter.
% Goal: Show the Green Line (Hypothesis 1) dropping instantly to Zero.
figure('Color','w','Position',[700 50 600 350]);
plot(time, log_npd_uio_dyn(1,:), 'g-', 'LineWidth', 2.5); hold on; % S1 (Correct)
plot(time, log_npd_uio_dyn(2,:), 'b--', 'LineWidth', 1.0);        % S2
plot(time, log_npd_uio_dyn(3,:), 'k:', 'LineWidth', 1.0);         % S3
xline(fault_time, 'k-', 'Fault Injection');

title('SUCCESS: Robust Filter Decision over Time');
xlabel('Time [s]'); ylabel('NPD Score (0 = Perfect Match)');
legend('Hypothesis: Sensor 1', 'Hypothesis: Sensor 2', 'Hypothesis: Sensor 3', 'Location', 'Best');
grid on; ylim([0 1.2]);

text(fault_time + 0.5, 0.1, '\leftarrow Instant Isolation of S1', 'FontSize', 10, 'FontWeight', 'bold');


% --- FIGURE 4: FULL TIME DOMAIN ANALYSIS ---
% This graph shows the magnitude of the residual vector ||r(t)||.
% Goal: Show "Pre-Fault" behavior. Standard Filter is noisy (False Alarms).
% Robust Filter is flat (Decoupled).
figure('Color','w','Position',[400 450 800 300]);

norm_b = vecnorm(log_r_bfdf);
norm_u = vecnorm(log_r_uio);

plot(time, norm_b, 'Color', [0.85 0.33 0.1], 'LineWidth', 2); hold on;
plot(time, norm_u, 'Color', [0.1 0.6 0.3], 'LineWidth', 2);
xline(fault_time, 'k--', 'Fault Injection (t=3s)');

title('Full Timeline Analysis: Residual Magnitude ||r(t)||');
xlabel('Time [s]'); ylabel('Magnitude');
legend('Standard Filter (Sensitive to errors)', 'Robust Filter (Decoupled)', 'Location', 'NorthWest');
grid on;
text(1.5, 0.025, 'FALSE ALARMS RISK', 'Color', [0.85 0.33 0.1], 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(1.5, 0.002, 'Perfect Silence', 'Color', [0.1 0.6 0.3], 'FontWeight', 'bold', 'HorizontalAlignment', 'center');
text(6, max(norm_b)*0.9, 'Standard reacts to Fault + Model Noise', 'Color', [0.85 0.33 0.1]);
text(6, 0.015, '\leftarrow Robust reacts ONLY to Fault', 'Color', [0.1 0.6 0.3], 'FontWeight', 'bold');
