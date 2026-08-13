% =========================================================
%  MASS SENSITIVITY ANALYSIS
%  ---------------------------------------------------------
%  Author: FARID ALIYEV
%  Date:   August, 2026
%  ---------------------------------------------------------
%  Dissertation: Fatty Acid PCM Cooling of PV Panels
%  Site: Jabrayil, Azerbaijan (39.33°N, 47.05°E)
%
%  Sweeps the PCM mass per panel (4 to 40 kg) at two melting
%  temperatures, 40°C (sweep optimum) and 44.13°C (the
%  selected MA-SA eutectic), to establish how the net energy
%  result depends on latent storage capacity. Model and all
%  other parameters as in Section 3.6 of the report.
%
%  HOW TO RUN:
%    1. Place this script and PVGIS_TMY_Summer_JJA.mat in
%       the same working directory (MATLAB R2020a or later).
%    2. Run:  Code3_Mass_Sensitivity
%       (14 cases; runtime a few minutes on a desktop PC)
%
%  INPUTS:  PVGIS_TMY_Summer_JJA.mat
%  OUTPUTS: mass_sensitivity_results.csv
% =========================================================

clc; clear; close all;

%% ── 1. LOAD DATA ─────────────────────────────────────────
DATA_PATH = 'PVGIS_TMY_Summer_JJA.mat';
fprintf('Loading PVGIS TMY summer data...\n');
d        = load(DATA_PATH);
G        = double(d.GHI(1,:));
T_amb    = double(d.T_amb(1,:));
v_wind   = double(d.wind_speed(1,:));
N        = length(G);
fprintf('  %d hourly records loaded.\n\n', N);

%% ── 2. FIXED PARAMETERS ─────────────────────────────────
A_pv    = 1.638;   m_pv  = 18.0;  Cp_pv = 900;
alpha   = 0.90;    eps_pv= 0.91;  eta_ref = 0.156;
beta_T  = 0.0045;  T_ref = 25.0;  sigma = 5.67e-8;

L_pcm   = 182400;  Cp_s  = 2090;  Cp_l  = 2200;
dT_melt = 4.0;
h_s     = 17;      h_m   = 200;   h_l   = 60;

dt      = 10;
n_sub   = 3600 / dt;

%% ── 3. BASELINE (runs once) ─────────────────────────────
fprintf('Running BASELINE simulation (no PCM)...\n');
P_B    = zeros(1, N);
T_pv_B = zeros(1, N);
T_pv_B_now = T_amb(1);
for i = 1:N
    G_i = G(i); Ta = T_amb(i); v = v_wind(i);
    h_f = 5.7 + 3.8*v;  h_b = 0.5*(5.7 + 3.8*v);
    T_sky_K = 0.0552*(Ta+273.15)^1.5;
    for s = 1:n_sub
        T_K = T_pv_B_now + 273.15;
        eta = max(eta_ref*(1-beta_T*(T_pv_B_now-T_ref)), 0);
        dT = (alpha*G_i*A_pv - eta*G_i*A_pv - h_f*A_pv*(T_pv_B_now-Ta) ...
              - eps_pv*sigma*A_pv*(T_K^4-T_sky_K^4) - h_b*A_pv*(T_pv_B_now-Ta)) ...
             /(m_pv*Cp_pv);
        T_pv_B_now = T_pv_B_now + dT*dt;
    end
    T_pv_B(i) = T_pv_B_now;
    P_B(i) = max(eta_ref*(1-beta_T*(T_pv_B_now-T_ref)),0)*G_i*A_pv;
end
E_baseline = sum(P_B)/1000;
idx_gen = G > 10;
fprintf('  Baseline summer energy: %.4f kWh/panel\n\n', E_baseline);

%% ── 4. MASS SWEEP AT TWO MELTING POINTS ─────────────────
m_range      = [4 8 12 16 20 30 40];          % kg
Tm_range     = [40.0, 44.13];                 % °C
n_m  = length(m_range);
n_tm = length(Tm_range);

results = zeros(n_m*n_tm, 8);  % Tm, m, E, gain, gain%, dT, %>=sol, minTpcm
row = 0;

fprintf('Starting mass sweep (%d cases)...\n\n', n_m*n_tm);

for k = 1:n_tm
    T_melt = Tm_range(k);
    T_sol  = T_melt - dT_melt/2;
    T_liq  = T_melt + dT_melt/2;

    for c = 1:n_m
        m_pcm = m_range(c);
        fprintf('[Tm=%.2f°C, m=%2d kg] ... ', T_melt, m_pcm);
        tic;

        T_pv_P_now = T_amb(1);
        T_pcm_now  = T_amb(1);
        H_pcm      = pcm_T_to_H(T_pcm_now, Cp_s, Cp_l, L_pcm, T_sol, T_liq);
        P_P    = zeros(1,N);
        T_pv_P = zeros(1,N);
        T_pcm_A= zeros(1,N);
        n_above= 0;

        for i = 1:N
            G_i = G(i); Ta = T_amb(i); v = v_wind(i);
            h_f = 5.7 + 3.8*v;  h_b = 0.5*(5.7 + 3.8*v);
            T_sky_K = 0.0552*(Ta+273.15)^1.5;
            for s = 1:n_sub
                T_K  = T_pv_P_now + 273.15;
                eta_p= max(eta_ref*(1-beta_T*(T_pv_P_now-T_ref)),0);
                if T_pcm_now < T_sol
                    h_eff = h_s;
                elseif T_pcm_now <= T_liq
                    h_eff = h_m;
                else
                    h_eff = h_l;
                end
                Q_to_pcm = h_eff*A_pv*(T_pv_P_now-T_pcm_now);
                dT_pv = (alpha*G_i*A_pv - eta_p*G_i*A_pv - h_f*A_pv*(T_pv_P_now-Ta) ...
                         - eps_pv*sigma*A_pv*(T_K^4-T_sky_K^4) - Q_to_pcm)/(m_pv*Cp_pv);
                T_pv_P_now = T_pv_P_now + dT_pv*dt;

                Q_pcm_loss = h_b*A_pv*(T_pcm_now-Ta);
                H_pcm = H_pcm + (Q_to_pcm - Q_pcm_loss)/m_pcm*dt;
                T_pcm_now = pcm_H_to_T(H_pcm,Cp_s,Cp_l,L_pcm,T_sol,T_liq);
            end
            T_pv_P(i)  = T_pv_P_now;
            T_pcm_A(i) = T_pcm_now;
            P_P(i)     = max(eta_ref*(1-beta_T*(T_pv_P_now-T_ref)),0)*G_i*A_pv;
            if T_pcm_now >= T_sol, n_above = n_above + 1; end
        end

        E_P   = sum(P_P)/1000;
        gain  = E_P - E_baseline;
        gpct  = gain/E_baseline*100;
        mdT   = mean(T_pv_B(idx_gen) - T_pv_P(idx_gen));
        row = row + 1;
        results(row,:) = [T_melt, m_pcm, E_P, gain, gpct, mdT, n_above/N*100, min(T_pcm_A)];

        fprintf('Done (%.1fs) | Gain: %+.4f kWh (%+.4f%%) | dT: %+.2f°C\n', ...
                toc, gain, gpct, mdT);
    end
    fprintf('\n');
end

%% ── 5. SUMMARY + CSV ────────────────────────────────────
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('  MASS SENSITIVITY RESULTS — baseline %.4f kWh\n', E_baseline);
fprintf('════════════════════════════════════════════════════════════════\n');
fprintf('  %-9s %-8s %-9s %-11s %-10s %-9s %-9s %-9s\n', ...
        'Tm(°C)','m(kg)','kg/m²','Gain(kWh)','Gain(%)','dT(°C)','>=sol(%)','Tpcm_min');
for r = 1:row
    fprintf('  %-9.2f %-8.0f %-9.1f %-11.4f %-10.4f %-9.2f %-9.1f %-9.1f\n', ...
            results(r,1), results(r,2), results(r,2)/A_pv, results(r,4), ...
            results(r,5), results(r,6), results(r,7), results(r,8));
end

fid = fopen('mass_sensitivity_results.csv','w');
fprintf(fid,'T_melt_C,m_pcm_kg,kg_per_m2,Baseline_kWh,PCM_kWh,Gain_kWh,Gain_pct,Mean_daytime_dT_C,Hrs_above_solidus_pct,Min_T_PCM_C\n');
for r = 1:row
    fprintf(fid,'%.2f,%.0f,%.2f,%.4f,%.4f,%.4f,%.4f,%.3f,%.2f,%.2f\n', ...
            results(r,1),results(r,2),results(r,2)/A_pv,E_baseline, ...
            results(r,3),results(r,4),results(r,5),results(r,6), ...
            results(r,7),results(r,8));
end
fclose(fid);
fprintf('\nResults saved to: mass_sensitivity_results.csv\n');

%% ── 6. FIGURE ───────────────────────────────────────────
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',11);
fig = figure('Position',[100 100 900 420],'Color','w');
hold on;
idx1 = results(1:n_m,:);
idx2 = results(n_m+1:2*n_m,:);
plot(m_range, idx1(:,5), 'o-', 'Color',[0.12 0.47 0.71],'LineWidth',1.8, ...
     'MarkerFaceColor',[0.12 0.47 0.71],'DisplayName','T_{melt} = 40°C (sweep optimum)');
plot(m_range, idx2(:,5), 's-', 'Color',[0.85 0.15 0.15],'LineWidth',1.8, ...
     'MarkerFaceColor',[0.85 0.15 0.15],'DisplayName','T_{melt} = 44.13°C (MA-SA eutectic)');
yline(0,'k-','LineWidth',1.0,'HandleVisibility','off');
xlabel('PCM Mass per Panel (kg)','FontSize',12);
ylabel('Summer Energy Gain (%)','FontSize',12);
title({'PCM Mass Sensitivity Sweep (4–40 kg per panel)';...
       'Jabrayil, Azerbaijan | JJA TMY | state-dependent h_{pvcm} (17/200/60 W/m²·K)'},...
      'FontSize',11,'FontWeight','normal');
legend('Location','best','FontSize',10);
grid on; box on;
xticks(m_range);
fprintf('Mass sensitivity figure generated.\n');

%% ── LOCAL FUNCTIONS ──────────────────────────────────────
function H = pcm_T_to_H(T, Cp_s, Cp_l, L, T_sol, T_liq)
    if T < T_sol
        H = Cp_s * T;
    elseif T <= T_liq
        H = Cp_s*T_sol + L*(T-T_sol)/(T_liq-T_sol);
    else
        H = Cp_s*T_sol + L + Cp_l*(T-T_liq);
    end
end

function T = pcm_H_to_T(H, Cp_s, Cp_l, L, T_sol, T_liq)
    H_sol = Cp_s*T_sol;
    H_liq = H_sol + L;
    if H < H_sol
        T = H/Cp_s;
    elseif H <= H_liq
        T = T_sol + (H-H_sol)/L*(T_liq-T_sol);
    else
        T = T_liq + (H-H_liq)/Cp_l;
    end
end
