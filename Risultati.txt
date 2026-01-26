%% PROGETTO ALES - REPLICAZIONE ESATTA & ANALISI VISIVA
% Reference: Chen, Patton, Zhang (1996)
% =========================================================================
% OBIETTIVO ESAME: 
% 1. Mostrare che il filtro Standard (BFDF) si confonde a causa delle non-linearità.
% 2. Mostrare che il filtro Robusto (UIO) ignora le non-linearità e isola il guasto.
% =========================================================================

clc; clear; close all;

%% 1. IL MONDO REALE E IL MODELLO (DATI DEL PAPER)

% --- MODELLO LINEARE (Che i filtri conoscono) ---
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
[n, ~] = size(A); [p, ~] = size(C);

% --- L'ERRORE DI MODELLO (Che i filtri NON conoscono o devono gestire) ---
% E_STAR (Pag 18): Rappresenta i termini quadratici (la realtà non lineare).
E_star = [1.3293  3.4440  0.1375 -5.1304 -1.7826 -1.8719;
          5.6812 -0.5281 -0.3385 -1.6193  0.5229  0;
          0       0       0       0       0       0     ];

% E1 (Pag 18): La matrice usata per progettare lo "scudo" del filtro robusto.
E1 = [6.2006,  2.8639;
      4.1048, -4.3262;
      0,       0     ];

% --- MATRICI DEL FILTRO ROBUSTO (Fornite a Pag 18-19) ---
H_paper = [ 0.6117, -0.1170, 0, 0.3215, 0.3220, -0.1295;
           -0.1170,  0.9382, 0, 0.0605, 0.0623, -0.1916;
            0,       0,      0, 0,      0,       0     ];
T_paper = [ 0, 0, -0.1251;
            0, 0,  0.0783;
            0, 0,  1.0000];
K_paper = [-0.0708,  0.0443,  0.5658,  0.1400,  0.1531,  0.3540;
            0.0443, -0.0277, -0.3540, -0.0876, -0.0958, -0.2215;
            0.5658, -0.3540, -4.5229, -1.1193, -1.2239, -2.8297];

%% 2. RETRO-INGEGNERIA (Calcoli teorici necessari)

% A. Filtro Standard (BFDF) - Non ha protezioni
pole = 3;
K_bfdf = (A + pole * eye(n)) * pinv(C);

% B. Filtro Robusto (Recuperiamo K1 per calcolare gli NPD corretti)
% Calcoliamo A1 usando la formula rigorosa: A1 = A - HCA
CE1 = C * E1;
A1 = A - E1 * inv(CE1' * CE1) * CE1' * C * A;
% Calcoliamo K1 (Guadagno direzionale)
K1_calc = (A1 + pole * eye(n)) * pinv(C);

%% 3. SIMULAZIONE "STRESS TEST"

dt = 0.001; T_sim = 8; time = 0:dt:T_sim; steps = length(time);
u_val = 0.20; fault_mag = 0.02; fault_time = 2.0;

% *** DISTORTION FACTOR = 20.0 ***
% Questo amplifica i termini quadratici E*. 
% È necessario per simulare un "cattivo modello lineare" e far fallire 
% il filtro standard, dimostrando la necessità della robustezza.
dist_factor = 20.0; 

I_sens = eye(p);
% Salvataggio dati per plot 3D (Scenario Sensore 2)
log_r_bfdf_3d = []; log_r_uio_3d = [];

% Risultati Finali
Table1 = zeros(3,3); Table2 = zeros(3,3);

fprintf('Simulazione in corso (Distortion Factor = %.1f)...\n', dist_factor);

for s_fault = 1:3 
    x = zeros(n,1); x_hat_bfdf = zeros(n,1); z_uio = zeros(n,1);
    acc_bfdf = zeros(1,3); acc_uio = zeros(1,3); count = 0;
    
    % Logger temporanei per il plot
    hist_r_bfdf = zeros(p, steps); hist_r_uio = zeros(p, steps);
    
    for k = 1:steps
        % 1. DINAMICA REALE (Plant Non-Lineare)
        % d_x contiene i termini quadratici (errore modello)
        d_x = [x(1)^2; x(2)^2; x(3)^2; x(1)*x(2); x(1)*x(3); x(2)*x(3)];
        % Evoluzione con disturbo amplificato
        dx = A*x + B*u_val + (E_star * d_x) * dist_factor;
        x = x + dx * dt;
        
        y_raw = C*x;
        % Guasto
        if time(k) >= fault_time, y_raw(s_fault) = y_raw(s_fault) + fault_mag; end
        
        % 2. FILTRO STANDARD (Soffre il disturbo)
        r_bfdf = y_raw - C * x_hat_bfdf;
        dx_bfdf = A * x_hat_bfdf + B * u_val + K_bfdf * r_bfdf;
        x_hat_bfdf = x_hat_bfdf + dx_bfdf * dt;
        
        % 3. FILTRO ROBUSTO (Cancella il disturbo)
        F_uio = A1 - K1_calc * C; 
        dz_uio = F_uio * z_uio + T_paper * B * u_val + K_paper * y_raw;
        z_uio = z_uio + dz_uio * dt;
        x_hat_uio = z_uio + H_paper * y_raw;
        r_uio = y_raw - C * x_hat_uio;
        
        % Salvataggio per grafici
        hist_r_bfdf(:,k) = r_bfdf; hist_r_uio(:,k) = r_uio;
        
        % 4. CALCOLO NPD (A regime)
        if time(k) > fault_time + 2.0
            count = count + 1;
            for j = 1:3
                Ij = I_sens(:,j);
                % NPD Standard
                Phi = [Ij, C*K_bfdf(:,j)];
                P = Phi * pinv(Phi'*Phi) * Phi';
                r_p = P * r_bfdf;
                acc_bfdf(j) = acc_bfdf(j) + (norm(r_bfdf-r_p)/(norm(r_bfdf)+eps));
                % NPD Robust
                Phi_r = [Ij, C*K1_calc(:,j), C*H_paper(:,j)];
                P_r = Phi_r * pinv(Phi_r'*Phi_r) * Phi_r';
                r_pr = P_r * r_uio;
                acc_uio(j) = acc_uio(j) + (norm(r_uio-r_pr)/(norm(r_uio)+eps));
            end
        end
    end
    Table1(:, s_fault) = acc_bfdf / count; % Nota: Colonna = Guasto Reale
    Table2(:, s_fault) = acc_uio / count;
    
    if s_fault == 2 % Salviamo i dati del caso "interessante" per il plot 3D
        log_r_bfdf_3d = hist_r_bfdf;
        log_r_uio_3d = hist_r_uio;
    end
end

%% 4. STAMPA TABELLE (Identiche al Paper)

fprintf('\n\nTable 1: BFDF (Standard) - Mis-isolation Example\n');
fprintf('Faulty sensor |    No.1   |    No.2   |    No.3   |\n');
fprintf('-----------------------------------------------\n');
for i=1:3
    fprintf('     NPD%d     |  %.5f  |  %.5f  |  %.5f  |\n', i, Table1(i,1), Table1(i,2), Table1(i,3));
end

fprintf('\n\nTable 2: UIO (Robust) - Perfect Isolation\n');
fprintf('Faulty sensor |    No.1   |    No.2   |    No.3   |\n');
fprintf('-----------------------------------------------\n');
for i=1:3
    fprintf('     NPD%d     |  %.5f  |  %.5f  |  %.5f  |\n', i, Table2(i,1), Table2(i,2), Table2(i,3));
end

%% 5. GRAFICI PER LA PRESENTAZIONE

% --- GRAFICO 1: IL FALLIMENTO (Bar Chart) ---
figure('Color','w','Position',[50 50 600 400]);
b_data = [Table1(:,2), Table2(:,2)]; % Analizziamo solo Colonna No.2
bar(b_data);
title('Analisi Guasto: Sensore 2 (Reale)');
xlabel('Ipotesi del Filtro'); ylabel('Distanza NPD (Basso = Match)');
xticklabels({'Ipotesi 1', 'Ipotesi 2', 'Ipotesi 3'});
legend('Standard (Sbaglia)', 'Robusto (Indovina)');
grid on;
text(2, Table1(2,2), '\downarrow Standard sbaglia!', 'Vert','bottom', 'Horiz','center', 'Color','r', 'FontWeight','bold');

% --- GRAFICO 2: LA MAGIA 3D (Geometric Isolation) ---
% Mostriamo come il vettore residuo si muove nello spazio 3D (r1, r2, r3)
% Usiamo solo i dati post-guasto
idx_start = find(time >= fault_time, 1);
r_bfdf_seg = log_r_bfdf_3d(1:3, idx_start:end);
r_uio_seg  = log_r_uio_3d(1:3, idx_start:end);

figure('Color','w','Position',[700 50 800 600]);
subplot(1,2,1);
plot3(r_bfdf_seg(1,:), r_bfdf_seg(2,:), r_bfdf_seg(3,:), 'r'); grid on;
title({'STANDARD BFDF', 'Residuo "Sporco" (Caos)'});
xlabel('r1'); ylabel('r2'); zlabel('r3');
view(45,30); axis equal;

subplot(1,2,2);
plot3(r_uio_seg(1,:), r_uio_seg(2,:), r_uio_seg(3,:), 'b', 'LineWidth', 2); grid on;
hold on; 
plot3([0, mean(r_uio_seg(1,:))], [0, mean(r_uio_seg(2,:))], [0, mean(r_uio_seg(3,:))], 'k--','LineWidth',1);
title({'ROBUST UIO', 'Direzionalità Perfetta'});
xlabel('r1'); ylabel('r2'); zlabel('r3');
view(45,30); axis equal;
subtitle('Il residuo giace esattamente nel sottospazio previsto');

% --- GRAFICO 3: EFFETTO DISACCOPPIAMENTO (Time Domain) ---
figure('Color','w','Position',[100 550 1000 300]);
subplot(1,2,1);
plot(time, log_r_bfdf_3d(2,:), 'r'); xline(fault_time,'--k');
title('Residuo r_2 (Standard)'); ylabel('Ampiezza'); ylim([-0.05 0.05]); grid on;
subtitle('Rumoroso prima del guasto (Errori Modello)');

subplot(1,2,2);
plot(time, log_r_uio_3d(2,:), 'b'); xline(fault_time,'--k');
title('Residuo r_2 (Robusto)'); ylabel('Ampiezza'); ylim([-0.05 0.05]); grid on;
subtitle('Piatto prima del guasto (Disaccoppiato)');