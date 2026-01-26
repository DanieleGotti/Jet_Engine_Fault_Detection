%% DESIGN OF UNKNOWN INPUT OBSERVERS AND ROBUST FAULT DETECTION FILTERS
% Reproduction of the paper by J. Chen, R.J. Patton, H.Y. Zhang
% Project: ALES Course - University of Bergamo
% -------------------------------------------------------------------------
clc; clear; close all;

%% 1. SYSTEM DEFINITION (From Page 17 & 18)
% Linearized Model Matrices
A = [-1.5581,  0.6925,  0.3974;
      0.2619, -2.2228,  0.2238;
      0       0      -10    ];

B = [0; 0; 10];

% Output Matrix (6 outputs, 3 states)
C = [1       0       0;
     0       1       0;
     0       0       1;
     0.55107 0.13320 0.30603;
     0.55217 0.13526 0.32912;
    -0.25693 -0.23625 0.61299];

% Identified Matrix E* for Quadratic Terms (Page 18)
% This represents the "True" nonlinear behavior (Modeling Errors)
E_star = [1.3293  3.4440  0.1375 -5.1304 -1.7826 -1.8719;
          5.6812 -0.5281 -0.3385 -1.6193  0.5229  0;
          0       0       0       0       0       0     ];

% Matrix E1 (Decomposition of E*) used for ROBUST Filter Design
% The filter knows only this direction of disturbance
E1 = [6.2006,  2.8639;
      4.1048, -4.3262;
      0,       0     ];

[n, m_in] = size(B);
[p, ~] = size(C);

%% 2. FILTER DESIGN

% --- A. Standard Beard Fault Detection Filter (BFDF) ---
% Design: Place eigenvalues at -3
pole_target = -3;
% Formula: K = (A + sigma*I) * pinv(C) (Page 14 eq 25 simplified logic)
% Note: Since C is not square invertible, we use pseudo-inverse for design
% The paper says "rank(C)=3... gain determined as K=(3I+A)C*" (Page 17)
K_bfdf = (A - pole_target * eye(n)) * pinv(C);

% --- B. Robust Unknown Input Observer (UIO) ---
% Step 1: Compute H to decouple disturbances (Lemma 1, Eq 12)
% H = E1 * inv((C*E1)' * (C*E1)) * (C*E1)'
CE1 = C * E1;
H_uio = E1 * inv(CE1' * CE1) * CE1';

% Step 2: Compute T (The "Shield" matrix) (Eq 6)
T_uio = eye(n) - H_uio * C;

% Step 3: Compute System Matrix A1 for the transformed state (Eq 13)
A1 = A - H_uio * C * A; % Note: The paper writes A1 = T*A conceptually in derivation
% Let's stick to Eq 13 explicitly: A1 = A - E[...]... which is A - HCA
% Or simply A1 = T*A (since T = I - HC)
A1 = T_uio * A;

% Step 4: Design K1 for directionality (Eq 25 applied to A1)
% "All eigenvalues of A1 - K1C set to -3"
K1 = (A1 - pole_target * eye(n)) * pinv(C);

% Step 5: Compute Final Observer Matrices (Eq 7, 8, 4)
F_uio = A1 - K1 * C;
K2 = F_uio * H_uio;
K_uio = K1 + K2; % Total Gain for y(t)

%% 3. SIMULATION LOOP (Reproducing Tables 1 & 2)

% Simulation settings
dt = 0.001;
T_final = 5;
time = 0:dt:T_final;
steps = length(time);
u_val = 0.20; % Input 20%

% Storage for results
NPD_Table_BFDF = zeros(3,3); % Rows: Fault Scenarios, Cols: NPD1, NPD2, NPD3
NPD_Table_UIO  = zeros(3,3);

% Fault Parameters
fault_mag = 0.02; % 2% offset
fault_time = 1.0; % Fault occurs at 1s
fault_idx_start = find(time >= fault_time, 1);

% Signature Directions / Subspaces for Isolation (Definition 6)
% For NPD calculation we need the projection matrices Phi_j
% Phi_j = [Ij, C*k_j, C*h_j]
I_mat = eye(p); 

disp('Starting Simulations...');

for scenario = 1:3 % Loop over faulty sensor 1, 2, 3
    
    fprintf('Simulating Fault in Sensor %d...\n', scenario);
    
    % Initialize States
    x = zeros(n,1);         % True System State
    x_hat_bfdf = zeros(n,1);% BFDF State
    z_uio = zeros(n,1);     % UIO State z
    
    % Data loggers for plotting (only for Scenario 2 as example)
    if scenario == 2
        log_r_bfdf = zeros(p, steps);
        log_r_uio  = zeros(p, steps);
    end
    
    % Variables to accumulate NPD stats (post-transient)
    sum_NPD_bfdf = zeros(1,3);
    sum_NPD_uio  = zeros(1,3);
    count_stats = 0;
    
    for k = 1:steps
        t = time(k);
        
        % --- A. REAL PLANT DYNAMICS (Nonlinear proxy) ---
        % Quadratic terms d(x) (Page 18)
        d_x = [x(1)^2; x(2)^2; x(3)^2; x(1)*x(2); x(1)*x(3); x(2)*x(3)];
        
        % State Update (Euler integration)
        dx = A*x + B*u_val + E_star * d_x; % The "Real" system has modeling errors!
        x = x + dx * dt;
        
        % Output Measurement
        y_clean = C*x;
        
        % Inject Sensor Fault
        y_meas = y_clean;
        if k >= fault_idx_start
            % Add 2% fault to the specific sensor (scenario)
            % Assuming nominal value around equilibrium is normalized to ~1 for simulation scaling
            % or simply additive 0.02 as per paper implication of "small variations"
            y_meas(scenario) = y_meas(scenario) + fault_mag; 
        end
        
        % --- B. BFDF FILTER UPDATE ---
        dy_hat_bfdf = C * x_hat_bfdf;
        r_bfdf = y_meas - dy_hat_bfdf;
        
        dx_hat_bfdf = A * x_hat_bfdf + B * u_val + K_bfdf * r_bfdf;
        x_hat_bfdf = x_hat_bfdf + dx_hat_bfdf * dt;
        
        % --- C. ROBUST UIO UPDATE ---
        % Observer dynamics: dz = F*z + T*B*u + K*y
        dz_uio = F_uio * z_uio + T_uio * B * u_val + K_uio * y_meas;
        z_uio = z_uio + dz_uio * dt;
        
        % State Reconstruction: x_hat = z + H*y
        x_hat_uio = z_uio + H_uio * y_meas;
        r_uio = y_meas - C * x_hat_uio;
        
        % --- D. LOGGING & NPD CALCULATION ---
        if scenario == 2
            log_r_bfdf(:,k) = r_bfdf;
            log_r_uio(:,k) = r_uio;
        end
        
        % Calculate NPD only after fault and transient (e.g., > 1.5s)
        if t > fault_time + 0.5
            count_stats = count_stats + 1;
            
            for j = 1:3 % Calculate NPD w.r.t Sensor 1, 2, 3 signatures
                % 1. Create Subspace Basis Phi_j
                Ij = I_mat(:,j);
                
                % BFDF Subspace: Span{Ij, C*K(:,j)} (Approx for standard filter)
                Phi_bfdf = [Ij, C*K_bfdf(:,j)]; 
                
                % UIO Subspace: Span{Ij, C*K1(:,j), C*H(:,j)} (Def 6, Eq 27)
                Phi_uio = [Ij, C*K1(:,j), C*H_uio(:,j)];
                
                % 2. Projection (Eq 29)
                % P = Phi * inv(Phi'*Phi) * Phi'
                % Robust inverse using pinv
                Proj_bfdf = Phi_bfdf * pinv(Phi_bfdf'*Phi_bfdf) * Phi_bfdf';
                Proj_uio  = Phi_uio  * pinv(Phi_uio'*Phi_uio)  * Phi_uio';
                
                r_star_bfdf = Proj_bfdf * r_bfdf;
                r_star_uio  = Proj_uio * r_uio;
                
                % 3. NPD (Eq 30)
                norm_r_bfdf = norm(r_bfdf);
                norm_r_uio  = norm(r_uio);
                
                val_bfdf = norm(r_bfdf - r_star_bfdf) / (norm_r_bfdf + 1e-6);
                val_uio  = norm(r_uio - r_star_uio)  / (norm_r_uio + 1e-6);
                
                sum_NPD_bfdf(j) = sum_NPD_bfdf(j) + val_bfdf;
                sum_NPD_uio(j)  = sum_NPD_uio(j)  + val_uio;
            end
        end
    end
    
    % Average NPDs
    NPD_Table_BFDF(scenario, :) = sum_NPD_bfdf / count_stats;
    NPD_Table_UIO(scenario, :)  = sum_NPD_uio  / count_stats;
end

%% 4. DISPLAY RESULTS (Reproducing Tables)

fprintf('\n=== TABLE 1: BFDF (Standard) ===\n');
fprintf('Faulty Sens | NPD1     | NPD2     | NPD3\n');
for i=1:3
    fprintf('Sensor %d    | %.5f  | %.5f  | %.5f\n', i, NPD_Table_BFDF(i,1), NPD_Table_BFDF(i,2), NPD_Table_BFDF(i,3));
end

fprintf('\n=== TABLE 2: ROBUST UIO (Paper Method) ===\n');
fprintf('Faulty Sens | NPD1     | NPD2     | NPD3\n');
for i=1:3
    fprintf('Sensor %d    | %.5f  | %.5f  | %.5f\n', i, NPD_Table_UIO(i,1), NPD_Table_UIO(i,2), NPD_Table_UIO(i,3));
end

%% 5. VISUALIZATION ("WOW" Factor)

% Plot 1: Residuals Time History (Scenario 2)
figure('Name', 'Residual Analysis: Fault Sensor 2', 'Color', 'w', 'Position', [100 100 1000 500]);
subplot(2,1,1);
plot(time, log_r_bfdf(1:3,:), 'LineWidth', 1.5);
xline(fault_time, '--r', 'Fault Injection');
title('Standard BFDF Residuals (Note: Noisy due to Modeling Errors)', 'FontSize', 12);
legend('r_1', 'r_2', 'r_3'); grid on;
ylabel('Amplitude');

subplot(2,1,2);
plot(time, log_r_uio(1:3,:), 'LineWidth', 1.5);
xline(fault_time, '--r', 'Fault Injection');
title('Robust UIO Residuals (Note: Disturbance Decoupled)', 'FontSize', 12);
legend('r_1', 'r_2', 'r_3'); grid on;
ylabel('Amplitude'); xlabel('Time (s)');

% Plot 2: NPD Comparison Bar Chart
figure('Name', 'NPD Performance Comparison', 'Color', 'w');
b = bar([NPD_Table_BFDF(2,:); NPD_Table_UIO(2,:)]');
xticklabels({'NPD_1', 'NPD_2', 'NPD_3'});
legend('Standard BFDF', 'Robust UIO');
title('Isolation Performance for Fault in Sensor 2');
ylabel('NPD (Lower is Better)');
grid on;
text(1:3, NPD_Table_BFDF(2,:), num2str(NPD_Table_BFDF(2,:)', '%.2f'), 'vert','bottom','horiz','center'); 
text(1:3, NPD_Table_UIO(2,:), num2str(NPD_Table_UIO(2,:)', '%.2f'), 'vert','bottom','horiz','center'); 

% Plot 3: 3D Trajectory of Residual Vector (Geometric Interpretation)
% We visualize the first 3 components of the residual vector
figure('Name', '3D Residual Trajectory', 'Color', 'w');
plot3(log_r_uio(1,:), log_r_uio(2,:), log_r_uio(3,:), 'b', 'LineWidth', 1.5); hold on;
plot3(log_r_uio(1,1), log_r_uio(2,1), log_r_uio(3,1), 'go', 'MarkerFaceColor','g'); % Start
plot3(log_r_uio(1,end), log_r_uio(2,end), log_r_uio(3,end), 'ro', 'MarkerFaceColor','r'); % End
grid on; axis equal;
xlabel('Residual 1'); ylabel('Residual 2'); zlabel('Residual 3');
title('Geometric Direction of Robust Residuals (Sensor 2 Fault)');
view(45, 30);
