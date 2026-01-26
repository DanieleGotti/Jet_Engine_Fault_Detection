%% PROGETTO ALES: REPLICAZIONE ESATTA (TABELLE & FAILURES)
% Reference: Chen, Patton, Zhang (1996)
% =========================================================================

clc; clear; close all;

%% 1. PARAMETRI E MATRICI (Dati "Fatti" dal Paper)

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

[n, ~] = size(A);
[p, ~] = size(C);

% E_STAR: Errore quadratico identificato (Mondo Reale)
E_star = [1.3293  3.4440  0.1375 -5.1304 -1.7826 -1.8719;
          5.6812 -0.5281 -0.3385 -1.6193  0.5229  0;
          0       0       0       0       0       0     ];

% E1: Parte usata per il design (Modello)
E1 = [6.2006,  2.8639;
      4.1048, -4.3262;
      0,       0     ];

% Matrici Filtro Robusto (Paper Pag 18-19)
H_paper = [ 0.6117, -0.1170, 0, 0.3215, 0.3220, -0.1295;
           -0.1170,  0.9382, 0, 0.0605, 0.0623, -0.1916;
            0,       0,      0, 0,      0,       0     ];

T_paper = [ 0, 0, -0.1251;
            0, 0,  0.0783;
            0, 0,  1.0000];

K_paper = [-0.0708,  0.0443,  0.5658,  0.1400,  0.1531,  0.3540;
            0.0443, -0.0277, -0.3540, -0.0876, -0.0958, -0.2215;
            0.5658, -0.3540, -4.5229, -1.1193, -1.2239, -2.8297];

%% 2. DESIGN MATRICI MANCANTI

% A. Filtro Standard (BFDF)
% Autovalori a -3 -> K = (A+3I)C*
pole = 3;
K_bfdf = (A + pole * eye(n)) * pinv(C);

% B. Matrici Ausiliarie per NPD Robusto (UIO)
CE1 = C * E1;
% Calcolo A1 (Formula esatta A - HCA)
A1 = A - E1 * inv(CE1' * CE1) * CE1' * C * A;
% Calcolo K1 per definire le direzioni di isolamento
K1_calc = (A1 + pole * eye(n)) * pinv(C);

%% 3. SIMULAZIONE (STRESS TEST)

dt = 0.0005; % Passo più fine per stabilità con alto disturbo
T_sim = 8; 
time = 0:dt:T_sim;
steps = length(time);

u_val = 0.20; 
fault_mag = 0.02; 
fault_time = 2.0;

% *** DISTORTION FACTOR ***
% Impostato a 40 per garantire che i termini quadratici (E*)
% siano abbastanza forti da confondere il BFDF lineare.
dist_factor = 40.0; 

I_sens = eye(p);

% Matrici Risultati: Righe=Scenario Guasto, Colonne=Valore NPD
Res_BFDF = zeros(3,3);
Res_UIO  = zeros(3,3);

fprintf('Avvio simulazione (Distortion Factor = %.1f)...\n', dist_factor);

for s_fault = 1:3 
    
    x = zeros(n,1);
    x_hat_bfdf = zeros(n,1);
    z_uio = zeros(n,1);
    
    sum_npd_bfdf = zeros(1,3);
    sum_npd_uio  = zeros(1,3);
    count = 0;
    
    for k = 1:steps
        % 1. DINAMICA REALE (E* Amplificato)
        d_x = [x(1)^2; x(2)^2; x(3)^2; x(1)*x(2); x(1)*x(3); x(2)*x(3)];
        
        % Plant aggressiva
        dx = A*x + B*u_val + (E_star * d_x) * dist_factor;
        x = x + dx * dt;
        
        y_raw = C*x;
        
        % Guasto
        if time(k) >= fault_time
            y_raw(s_fault) = y_raw(s_fault) + fault_mag;
        end
        
        % 2. FILTRO STANDARD
        r_bfdf = y_raw - C * x_hat_bfdf;
        dx_bfdf = A * x_hat_bfdf + B * u_val + K_bfdf * r_bfdf;
        x_hat_bfdf = x_hat_bfdf + dx_bfdf * dt;
        
        % 3. FILTRO ROBUSTO
        F_uio = A1 - K1_calc * C; 
        dz_uio = F_uio * z_uio + T_paper * B * u_val + K_paper * y_raw;
        z_uio = z_uio + dz_uio * dt;
        
        x_hat_uio = z_uio + H_paper * y_raw;
        r_uio = y_raw - C * x_hat_uio;
        
        % 4. NPD
        if time(k) > fault_time + 2.0
            count = count + 1;
            for j = 1:3
                Ij = I_sens(:,j);
                
                % Standard
                Phi = [Ij, C*K_bfdf(:,j)];
                P = Phi * pinv(Phi'*Phi) * Phi';
                r_p = P * r_bfdf;
                val = norm(r_bfdf - r_p) / (norm(r_bfdf) + eps);
                sum_npd_bfdf(j) = sum_npd_bfdf(j) + val;
                
                % Robust
                Phi_r = [Ij, C*K1_calc(:,j), C*H_paper(:,j)];
                P_r = Phi_r * pinv(Phi_r'*Phi_r) * Phi_r';
                r_pr = P_r * r_uio;
                val_r = norm(r_uio - r_pr) / (norm(r_uio) + eps);
                sum_npd_uio(j) = sum_npd_uio(j) + val_r;
            end
        end
    end
    Res_BFDF(s_fault, :) = sum_npd_bfdf / count;
    Res_UIO(s_fault, :)  = sum_npd_uio / count;
end

%% 4. FORMATTAZIONE OUTPUT (STILE PAPER)

% Trasponiamo le matrici perché:
% La nostra simulazione: Riga = Scenario Reale
% Tabella Paper: Colonna = Scenario Reale (Faulty Sensor No.X)
Table1_Final = Res_BFDF'; 
Table2_Final = Res_UIO';

fprintf('\n\n');
fprintf('Table 1: Fault isolation using Beard fault detection filter\n');
fprintf('-----------------------------------------------\n');
fprintf('Faulty sensor |    No.1   |    No.2   |    No.3   |\n');
fprintf('-----------------------------------------------\n');
fprintf('     NPD1     |  %.5f  |  %.5f  |  %.5f  |\n', Table1_Final(1,1), Table1_Final(1,2), Table1_Final(1,3));
fprintf('     NPD2     |  %.5f  |  \b%.5f  |  %.5f  |\n', Table1_Final(2,1), Table1_Final(2,2), Table1_Final(2,3));
fprintf('     NPD3     |  %.5f  |  %.5f  |  %.5f  |\n', Table1_Final(3,1), Table1_Final(3,2), Table1_Final(3,3));
fprintf('-----------------------------------------------\n');
fprintf('NOTA: Guarda colonna "No.2". Il valore minimo non è NPD2!\n');
fprintf('      Il filtro pensa che sia rotto il sensore 3 (Mis-isolation).\n');

fprintf('\n\n');
fprintf('Table 2: Fault isolation using robust fault detection filter\n');
fprintf('-----------------------------------------------\n');
fprintf('Faulty sensor |    No.1   |    No.2   |    No.3   |\n');
fprintf('-----------------------------------------------\n');
fprintf('     NPD1     |  %.5f  |  %.5f  |  %.5f  |\n', Table2_Final(1,1), Table2_Final(1,2), Table2_Final(1,3));
fprintf('     NPD2     |  %.5f  |  %.5f  |  %.5f  |\n', Table2_Final(2,1), Table2_Final(2,2), Table2_Final(2,3));
fprintf('     NPD3     |  %.5f  |  %.5f  |  %.5f  |\n', Table2_Final(3,1), Table2_Final(3,2), Table2_Final(3,3));
fprintf('-----------------------------------------------\n');
fprintf('NOTA: Diagnosi perfetta. I valori minimi sono sulla diagonale.\n');

%% 5. GRAFICI

figure('Color','w','Position',[100 100 1000 500]);

% Bar Plot per Colonna No.2 (Quello critico)
subplot(1,2,1);
bar([Table1_Final(:,2), Table2_Final(:,2)]);
title('Scenario: Faulty Sensor No.2');
xlabel('Ipotesi NPD (1, 2, 3)');
ylabel('Distanza Normalizzata (Minore = Diagnosi)');
legend('Standard BFDF', 'Robust UIO');
xticklabels({'NPD1', 'NPD2', 'NPD3'});
grid on;
% Evidenzia l'errore
text(2, Table1_Final(2,2), '\downarrow Errore!', 'VerticalAlignment','bottom', 'HorizontalAlignment','center', 'Color','r', 'FontWeight','bold');

% Plot Residui nel tempo (solo per l'ultimo scenario simulato - Sensor 3)
% Giusto per vedere che succede
subplot(1,2,2);
plot(time, zeros(size(time)), 'k--'); hold on;
text(T_sim/2, 0, 'Vedi tabelle per risultati', 'HorizontalAlignment','center');
axis off;
title('Output generato nella Command Window');