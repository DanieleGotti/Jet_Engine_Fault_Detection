%% DESIGN OF UNKNOWN INPUT OBSERVERS - REPLICAZIONE ESATTA
% ALES Course Project
% Reference: Chen, Patton, Zhang (1996)
% -------------------------------------------------------------------------
clc; clear; close all;

%% 1. DEFINIZIONE DEL SISTEMA E DATI DAL PAPER

% --- Matrici del Modello Lineare (Pagina 17) ---
A = [-1.5581,  0.6925,  0.3974;
      0.2619, -2.2228,  0.2238;
      0       0      -10    ];

B = [0; 0; 10];

C = [1       0       0;
     0       1       0;
     0       0       1;
     0.55107 0.13320 0.30603;
     0.55217 0.13526 0.32912;
    -0.25693 -0.23625 0.61299];

[n, m_in] = size(B);
[p, ~] = size(C);

% --- Matrice di Errore del Modello E* (Pagina 18) ---
% Questa matrice mappa i termini quadratici sugli stati. 
% È stata identificata dagli autori via Least Squares.
E_star = [1.3293  3.4440  0.1375 -5.1304 -1.7826 -1.8719;
          5.6812 -0.5281 -0.3385 -1.6193  0.5229  0;
          0       0       0       0       0       0     ];

% --- Matrici del Filtro Robusto Fornite dal Paper (Pagina 18-19) ---
% NON le calcoliamo, usiamo quelle degli autori per replicare i risultati.
H_paper = [ 0.6117, -0.1170, 0, 0.3215, 0.3220, -0.1295;
           -0.1170,  0.9382, 0, 0.0605, 0.0623, -0.1916;
            0,       0,      0, 0,      0,       0     ];

T_paper = [ 0, 0, -0.1251;
            0, 0,  0.0783;
            0, 0,  1.0000]; % Nota: T è 3x3

% Matrice K del filtro robusto (Pagina 19)
K_paper = [-0.0708,  0.0443,  0.5658,  0.1400,  0.1531,  0.3540;
            0.0443, -0.0277, -0.3540, -0.0876, -0.0958, -0.2215;
            0.5658, -0.3540, -4.5229, -1.1193, -1.2239, -2.8297];

%% 2. CALCOLI MANCANTI (Quelli che dobbiamo fare noi)

% --- A. Matrice K per il Filtro Standard (BFDF - Tabella 1) ---
% Il paper dice a Pag 17: "K = (3I + A)C*" 
% Usiamo la pseudo-inversa (pinv) per C*
sigma = 3; 
K_bfdf = (A + sigma * eye(n)) * pinv(C);

% --- B. Ricostruzione di K1 per il calcolo NPD Robusto ---
% Ci serve K1 per definire le direzioni di isolamento (Eq. 27).
% Il paper non da K1 esplicito, ma dice come calcolarlo (Eq. 25).
% Prima calcoliamo A1 = T*A (Eq. 13 semplificata o derivata)
% Verifica: A1 dovrebbe essere stabile o stabilizzabile.
A1 = T_paper * A; % DA VERIFICARE
% Calcolo K1 (autovalori a -3 come da testo pag 18)
K1_calc = (A1 + sigma * eye(n)) * pinv(C);

% Verifica di coerenza (Opzionale): K_paper dovrebbe essere circa K1 + F*H
% Ma usiamo K_paper per la simulazione per fedeltà.

%% 3. SIMULAZIONE

dt = 0.001;
T_sim = 10; % Tempo simulazione esteso per steady state
time = 0:dt:T_sim;
steps = length(time);
u_val = 0.20; 

% Parametri Guasto
fault_mag = 0.02; 
fault_time = 2.0; % Guasto a 2s
idx_fault = find(time >= fault_time, 1);

% Tabelle risultati
NPD_Res_BFDF = zeros(3,3);
NPD_Res_UIO  = zeros(3,3);

% Matrice Identità per le direzioni di guasto sensori (Ij)
I_sens = eye(p);

fprintf('Avvio simulazione...\n');

for s_fault = 1:3 % Per ogni sensore guasto (1, 2, 3)
    
    % Inizializza stati
    x = zeros(n,1);
    x_hat_bfdf = zeros(n,1);
    z_uio = zeros(n,1);
    
    % Accumulatori per media NPD
    npd_acc_bfdf = zeros(1,3);
    npd_acc_uio = zeros(1,3);
    count = 0;
    
    % Logger per plot
    if s_fault == 2
        log_r_bfdf = zeros(p, steps);
        log_r_uio = zeros(p, steps);
    end

    for k = 1:steps
        % --- 1. DINAMICA REALE (Plant con termini quadratici) ---
        % Vettore d(x) (Pagina 18)
        d_x = [x(1)^2; x(2)^2; x(3)^2; x(1)*x(2); x(1)*x(3); x(2)*x(3)];
        
        % Evoluzione stato vero: Ax + Bu + E*d(x)
        dx = A*x + B*u_val + E_star * d_x;
        x = x + dx * dt;
        
        % Uscita Misurata
        y = C*x;
        
        % Iniezione Guasto
        if k >= idx_fault
            y(s_fault) = y(s_fault) + fault_mag;
        end
        
        % --- 2. FILTRO STANDARD (BFDF) ---
        % Observer classico: dx_hat = A x_hat + B u + K(y - C x_hat)
        y_hat_bfdf = C * x_hat_bfdf;
        r_bfdf = y - y_hat_bfdf;
        dx_hat_bfdf = A * x_hat_bfdf + B * u_val + K_bfdf * r_bfdf;
        x_hat_bfdf = x_hat_bfdf + dx_hat_bfdf * dt;
        
        % --- 3. FILTRO ROBUSTO (UIO - Paper) ---
        % Eq. 2: dz = Fz + TBu + Ky
        % Ma ci manca F esplicita nel paper. 
        % La calcoliamo da Eq 7: F = A - HCA - K1C = A1 - K1C
        F_uio = A1 - K1_calc * C; 
        
        % Nota: Usiamo K_paper fornita per il termine Ky, e T_paper per TBu
        dz_uio = F_uio * z_uio + T_paper * B * u_val + K_paper * y;
        z_uio = z_uio + dz_uio * dt;
        
        % Ricostruzione stato e residuo
        x_hat_uio = z_uio + H_paper * y;
        r_uio = y - C * x_hat_uio;
        
        % --- LOGGING & NPD ---
        if s_fault == 2
            log_r_bfdf(:,k) = r_bfdf;
            log_r_uio(:,k) = r_uio;
        end
        
        % Calcolo NPD (solo a regime post-guasto)
        if time(k) > fault_time + 2.0 
            count = count + 1;
            for j = 1:3 % Contro le 3 firme di guasto
                % Direzioni (Definizione 6 e dintorni)
                Ij = I_sens(:,j);
                
                % --- Subspace BFDF ---
                % Span{Ij, C*K_bfdf(:,j)} 
                % Nota: per BFDF standard spesso si usa solo la direzione C*Kj per attuatori,
                % ma per sensori il residuo giace nel piano definito da Ij e C*Kj
                Phi_bfdf = [Ij, C*K_bfdf(:,j)];
                Proj_bfdf = Phi_bfdf * pinv(Phi_bfdf' * Phi_bfdf) * Phi_bfdf';
                r_proj_bfdf = Proj_bfdf * r_bfdf;
                val_bfdf = norm(r_bfdf - r_proj_bfdf) / (norm(r_bfdf) + eps);
                npd_acc_bfdf(j) = npd_acc_bfdf(j) + val_bfdf;
                
                % --- Subspace UIO ---
                % Span{Ij, C*K1(:,j), C*H(:,j)} (Eq 27 & Fig 2)
                Phi_uio = [Ij, C*K1_calc(:,j), C*H_paper(:,j)];
                Proj_uio = Phi_uio * pinv(Phi_uio' * Phi_uio) * Phi_uio';
                r_proj_uio = Proj_uio * r_uio;
                val_uio = norm(r_uio - r_proj_uio) / (norm(r_uio) + eps);
                npd_acc_uio(j) = npd_acc_uio(j) + val_uio;
            end
        end
    end
    
    NPD_Res_BFDF(s_fault, :) = npd_acc_bfdf / count;
    NPD_Res_UIO(s_fault, :)  = npd_acc_uio / count;
end

%% 4. RISULTATI FINALI E PLOT

fprintf('\n\n--- TABELLA 1 (BFDF Standard) ---\n');
disp(array2table(NPD_Res_BFDF, 'VariableNames', {'NPD1','NPD2','NPD3'}, 'RowNames', {'Fault1','Fault2','Fault3'}));

fprintf('\n--- TABELLA 2 (UIO Robusto) ---\n');
disp(array2table(NPD_Res_UIO, 'VariableNames', {'NPD1','NPD2','NPD3'}, 'RowNames', {'Fault1','Fault2','Fault3'}));

% --- Grafici Extra ---
figure('Color','w','Position',[100 100 1000 600]);

% 1. Residui Standard
subplot(2,2,1);
plot(time, log_r_bfdf(1:3,:), 'LineWidth',1); 
xline(fault_time,'--k');
title('Residui BFDF (Guasto Sensore 2)');
xlabel('Tempo (s)'); grid on; legend('r1','r2','r3');

% 2. Residui Robusti
subplot(2,2,3);
plot(time, log_r_uio(1:3,:), 'LineWidth',1);
xline(fault_time,'--k');
title('Residui UIO Robusti (Guasto Sensore 2)');
xlabel('Tempo (s)'); grid on; legend('r1','r2','r3');

% 3. Confronto NPD a barre
subplot(2,2,[2,4]);
bar_data = [NPD_Res_BFDF(2,:); NPD_Res_UIO(2,:)]';
b = bar(bar_data);
xticklabels({'NPD1 (Wrong)','NPD2 (Correct)','NPD3 (Wrong)'});
legend('Standard','Robust');
title('Confronto Isolamento Guasto Sensore 2');
grid on;
ylabel('Normalized Projection Distance');

% 3D Trajectory
figure('Color','w');
plot3(log_r_uio(1,:), log_r_uio(2,:), log_r_uio(3,:), 'b'); hold on;
plot3(log_r_uio(1,end), log_r_uio(2,end), log_r_uio(3,end), 'ro', 'MarkerSize',10,'MarkerFaceColor','r');
grid on; xlabel('r1'); ylabel('r2'); zlabel('r3');
title('Traiettoria Residuo Robusto nello Spazio 3D');
view(45,30);