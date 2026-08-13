% =========================================================
%  MELTING-POINT SENSITIVITY ANALYSIS
%  ---------------------------------------------------------
%  Author: FARID ALIYEV
%  Date:   August, 2026
%  ---------------------------------------------------------
%  Dissertation: Fatty Acid PCM Cooling of PV Panels
%  Site: Jabrayil, Azerbaijan (39.33°N, 47.05°E)
%
%  Sweeps the PCM melting temperature from 30°C to 62°C in
%  2°C steps at a fixed PCM mass of 4 kg per panel, using
%  the heat transfer model of Section 3.6 of the report.
%  The lower part of the range (30–42°C) covers the
%  lower-melting eutectics reviewed in Section 2.5
%  (e.g. LA-MA 34.2°C, LA-SA 37°C).
%
%  NOTE: latent heat and specific heats are held at the
%  MA-SA 64:36 values for all melting points, so the sweep
%  isolates the effect of the melting temperature alone.
%
%  HOW TO RUN:
%    1. Place this script and PVGIS_TMY_Summer_JJA.mat in
%       the same working directory (MATLAB R2020a or later).
%    2. Run:  Code2_MeltingPoint_Sensitivity
%       (17 cases; runtime a few minutes on a desktop PC)
%
%  INPUTS:  PVGIS_TMY_Summer_JJA.mat
%  OUTPUTS: sensitivity_results.csv
%
%  FIXED PARAMETERS:
%    h_pvcm: state-dependent (solid 17 / melting 200 /
%            liquid 60 W/m²·K)
%    m_pcm  = 4.0 kg
%    dt     = 10 s
% =========================================================

clc; clear; close all;

%% ── 1. DATA PATH ─────────────────────────────────────────
DATA_PATH = 'PVGIS_TMY_Summer_JJA.mat';

%% ── 2. LOAD PVGIS CLIMATE DATA ───────────────────────────
fprintf('Loading PVGIS TMY summer data...\n');
d        = load(DATA_PATH);
G        = double(d.GHI(1,:));
T_amb    = double(d.T_amb(1,:));
v_wind   = double(d.wind_speed(1,:));
month    = double(d.month(1,:));
N        = length(G);
fprintf('  %d hourly records loaded.\n\n', N);

%% ── 3. FIXED PV PANEL PARAMETERS ────────────────────────
A_pv    = 1.638;
m_pv    = 18.0;
Cp_pv   = 900;
alpha   = 0.90;
eps_pv  = 0.91;
eta_ref = 0.156;
beta_T  = 0.0045;
T_ref   = 25.0;
sigma   = 5.67e-8;

%% ── 4. FIXED PCM PARAMETERS ─────────────────────────────
L_pcm   = 182400;     % Latent heat, MA-SA 64:36 [J/kg], held fixed for all cases
Cp_s    = 2090;       % Solid specific heat [J/kg·K]
Cp_l    = 2200;       % Liquid specific heat [J/kg·K]
m_pcm   = 4.0;        % PCM mass [kg]
dT_melt = 4.0;        % Mushy zone width [°C]
h_s     = 17;         % Solid phase conductance [W/m²·K]  — k_pcm/L = 0.170/0.010
h_m     = 200;        % Melting phase conductance [W/m²·K] — natural convection enhanced
h_l     = 60;         % Liquid phase conductance [W/m²·K]  — moderate convection
                      % (Huang et al., 2004; Hasan et al., 2014)

%% ── 5. TIME STEPPING ─────────────────────────────────────
dt      = 10;
n_sub   = 3600 / dt;  % 360 sub-steps per hour

%% ── 6. BASELINE SIMULATION (runs once) ──────────────────
fprintf('Running BASELINE simulation (no PCM)...\n');
P_B    = zeros(1, N);
T_pv_B = zeros(1, N);        % Baseline PV temperature (stored for ΔT metric)
T_pv_B_now = T_amb(1);

for i = 1:N
    G_i = G(i); Ta = T_amb(i); v = v_wind(i);
    h_f = 5.7 + 3.8*v;
    h_b = 0.5*(5.7 + 3.8*v);
    T_sky_K = 0.0552*(Ta+273.15)^1.5;
    for s = 1:n_sub
        T_K = T_pv_B_now + 273.15;
        eta = max(eta_ref*(1-beta_T*(T_pv_B_now-T_ref)), 0);
        Q_sol   = alpha*G_i*A_pv;
        Q_elec  = eta*G_i*A_pv;
        Q_conv_f= h_f*A_pv*(T_pv_B_now-Ta);
        Q_rad   = eps_pv*sigma*A_pv*(T_K^4-T_sky_K^4);
        Q_conv_b= h_b*A_pv*(T_pv_B_now-Ta);
        dT = (Q_sol-Q_elec-Q_conv_f-Q_rad-Q_conv_b)/(m_pv*Cp_pv);
        T_pv_B_now = T_pv_B_now + dT*dt;
    end
    T_pv_B(i) = T_pv_B_now;
    P_B(i) = max(eta_ref*(1-beta_T*(T_pv_B_now-T_ref)),0)*G_i*A_pv;
end
E_baseline = sum(P_B)/1000;  % [kWh]
idx_gen = G > 10;            % Daytime (generating) hours mask
fprintf('  Baseline summer energy: %.3f kWh/panel\n', E_baseline);
fprintf('  Baseline mean daytime panel temp: %.2f °C\n\n', mean(T_pv_B(idx_gen)));

%% ── 7. SENSITIVITY SWEEP ─────────────────────────────────
T_melt_range = 30:2:62;   % °C — 17 values, spanning the candidate eutectics
n_cases      = length(T_melt_range);

E_PCM        = zeros(1, n_cases);
E_gain       = zeros(1, n_cases);
E_gain_pct   = zeros(1, n_cases);
peak_T_PCM   = zeros(1, n_cases);
min_T_PCM    = zeros(1, n_cases);   % Min PCM temp — checks nightly re-solidification
pct_above_sol= zeros(1, n_cases);   % % of hours with T_pcm ≥ solidus (latent zone or liquid)
pct_mushy    = zeros(1, n_cases);   % % of hours inside the mushy zone (actively melting)
mean_dT      = zeros(1, n_cases);   % Mean daytime panel temp reduction vs baseline [°C]

fprintf('Starting sensitivity sweep (%d melting points, %d–%d°C)...\n', ...
        n_cases, T_melt_range(1), T_melt_range(end));
fprintf('\n');

for c = 1:n_cases
    T_melt = T_melt_range(c);
    T_sol  = T_melt - dT_melt/2;
    T_liq  = T_melt + dT_melt/2;

    fprintf('[%2d/%d] T_melt = %.0f°C (T_sol=%.1f, T_liq=%.1f)... ', ...
            c, n_cases, T_melt, T_sol, T_liq);
    tic;

    % Initialise state
    T_pv_P_now = T_amb(1);
    T_pcm_now  = T_amb(1);
    H_pcm      = pcm_T_to_H(T_pcm_now, Cp_s, Cp_l, L_pcm, T_sol, T_liq);

    P_P     = zeros(1,N);
    T_pv_P  = zeros(1,N);
    T_pcm_A = zeros(1,N);
    n_above = 0;
    n_mushy = 0;

    for i = 1:N
        G_i = G(i); Ta = T_amb(i); v = v_wind(i);
        h_f = 5.7 + 3.8*v;
        h_b = 0.5*(5.7 + 3.8*v);
        T_sky_K = 0.0552*(Ta+273.15)^1.5;

        for s = 1:n_sub
            T_K  = T_pv_P_now + 273.15;
            eta_p= max(eta_ref*(1-beta_T*(T_pv_P_now-T_ref)),0);
            Q_sol    = alpha*G_i*A_pv;
            Q_elec_p = eta_p*G_i*A_pv;
            Q_conv_f = h_f*A_pv*(T_pv_P_now-Ta);
            Q_rad    = eps_pv*sigma*A_pv*(T_K^4-T_sky_K^4);
            if T_pcm_now < T_sol
                h_eff = h_s;
            elseif T_pcm_now <= T_liq
                h_eff = h_m;
            else
                h_eff = h_l;
            end
            Q_to_pcm = h_eff*A_pv*(T_pv_P_now-T_pcm_now);

            dT_pv = (Q_sol-Q_elec_p-Q_conv_f-Q_rad-Q_to_pcm)/(m_pv*Cp_pv);
            T_pv_P_now = T_pv_P_now + dT_pv*dt;

            Q_pcm_loss = h_b*A_pv*(T_pcm_now-Ta);
            dH = (Q_to_pcm - Q_pcm_loss)/m_pcm;
            H_pcm = H_pcm + dH*dt;
            T_pcm_now = pcm_H_to_T(H_pcm,Cp_s,Cp_l,L_pcm,T_sol,T_liq);
        end  % end sub-step loop

        T_pv_P(i)  = T_pv_P_now;
        T_pcm_A(i) = T_pcm_now;
        P_P(i)     = max(eta_ref*(1-beta_T*(T_pv_P_now-T_ref)),0)*G_i*A_pv;
        if T_pcm_now >= T_sol, n_above = n_above + 1; end
        if T_pcm_now >= T_sol && T_pcm_now <= T_liq, n_mushy = n_mushy + 1; end
    end

    E_PCM(c)       = sum(P_P)/1000;
    E_gain(c)      = E_PCM(c) - E_baseline;
    E_gain_pct(c)  = (E_gain(c)/E_baseline)*100;
    peak_T_PCM(c)  = max(T_pcm_A);
    min_T_PCM(c)   = min(T_pcm_A);
    pct_above_sol(c) = (n_above/N)*100;
    pct_mushy(c)   = (n_mushy/N)*100;
    mean_dT(c)     = mean(T_pv_B(idx_gen) - T_pv_P(idx_gen));  % + = PCM cools panel

    elapsed = toc;
    fprintf('Done (%.1fs) | Gain: %+.4f kWh (%+.4f%%) | Mean daytime dT: %+.2f°C | >=solidus: %.1f%% of hrs\n', ...
            elapsed, E_gain(c), E_gain_pct(c), mean_dT(c), pct_above_sol(c));
end

%% ── 8. PRINT SUMMARY TABLE ───────────────────────────────
fprintf('\n');
fprintf('════════════════════════════════════════════════════════════════════════════════\n');
fprintf('  SENSITIVITY ANALYSIS RESULTS — Jabrayil, Azerbaijan (JJA TMY)\n');
fprintf('  Fixed: state-dependent h_pvcm (s=%.0f/m=%.0f/l=%.0f W/m²·K) | m_pcm=%.1f kg | dt=%.0fs\n', ...
        h_s, h_m, h_l, m_pcm, dt);
fprintf('════════════════════════════════════════════════════════════════════════════════\n');
fprintf('  %-9s %-11s %-11s %-11s %-9s %-9s %-9s %-9s %-9s\n', ...
        'Tmelt(°C)','PCM(kWh)','Gain(kWh)','Gain(%)','dT(°C)','>=sol(%)','mushy(%)','Tpcm_max','Tpcm_min');
fprintf('  %s\n', repmat('-',1,92));
for c = 1:n_cases
    fprintf('  %-9.0f %-11.4f %-11.4f %-11.4f %-9.2f %-9.1f %-9.1f %-9.1f %-9.1f\n', ...
            T_melt_range(c), E_PCM(c), E_gain(c), E_gain_pct(c), ...
            mean_dT(c), pct_above_sol(c), pct_mushy(c), peak_T_PCM(c), min_T_PCM(c));
end
fprintf('  %s\n', repmat('-',1,92));
fprintf('  Baseline energy: %.4f kWh | Baseline mean daytime T: %.2f°C\n', ...
        E_baseline, mean(T_pv_B(idx_gen)));
[best_gain, best_idx] = max(E_gain_pct);
fprintf('  LEAST-UNFAVOURABLE T_melt: %.0f°C  (gain = %+.4f%%)\n', ...
        T_melt_range(best_idx), best_gain);
fprintf('  NOTE: If Tpcm_min > T_sol for a case, the PCM is not fully re-solidifying\n');
fprintf('        overnight — latent capacity is partially unavailable the next day.\n\n');

%% ── 9. EXPORT TO CSV ────────────────────────────────────
fid = fopen('sensitivity_results.csv','w');
fprintf(fid,'T_melt_C,Baseline_kWh,PCM_kWh,Gain_kWh,Gain_pct,Mean_daytime_dT_C,Hrs_above_solidus_pct,Hrs_mushy_pct,Peak_T_PCM_C,Min_T_PCM_C\n');
for c = 1:n_cases
    fprintf(fid,'%.1f,%.4f,%.4f,%.4f,%.4f,%.3f,%.2f,%.2f,%.2f,%.2f\n', ...
            T_melt_range(c),E_baseline,E_PCM(c),E_gain(c), ...
            E_gain_pct(c),mean_dT(c),pct_above_sol(c),pct_mushy(c), ...
            peak_T_PCM(c),min_T_PCM(c));
end
fclose(fid);
fprintf('Results saved to: sensitivity_results.csv\n\n');

%% ── 10. PUBLICATION FIGURE ──────────────────────────────
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',11);

fig = figure('Position',[100 100 980 440],'Color','w');

% Left axis — energy gain %
yyaxis left
bar(T_melt_range, E_gain_pct, 0.5, 'FaceColor',[0.12 0.47 0.71],'FaceAlpha',0.8);
hold on;
plot(T_melt_range(best_idx), E_gain_pct(best_idx), 'r*', 'MarkerSize',12,'LineWidth',2);
yline(0,'k-','LineWidth',1.0,'HandleVisibility','off');
xlabel('PCM Melting Temperature (°C)','FontSize',12);
ylabel('Summer Energy Gain (%)','FontSize',12);
ylim_left = max(abs(E_gain_pct))*1.4;
ylim([-ylim_left ylim_left]);

% Candidate eutectic markers (Section 2.5 table)
cand_T   = [34.2, 37.0, 44.13];
cand_lbl = {'LA-MA (34.2°C)','LA-SA (37°C)','MA-SA (44.1°C)'};
for k = 1:3
    xline(cand_T(k), ':', 'Color',[0.3 0.3 0.3], 'LineWidth',1.1, 'HandleVisibility','off');
    text(cand_T(k), -ylim_left*0.85, cand_lbl{k}, 'Rotation',90, ...
         'FontSize',8,'FontName','Times New Roman','Color',[0.3 0.3 0.3], ...
         'HorizontalAlignment','left');
end

% Right axis — % hours at/above solidus
yyaxis right
plot(T_melt_range, pct_above_sol, 'ko-','LineWidth',1.5,'MarkerFaceColor','k','MarkerSize',5);
ylabel('Hours with PCM at/above Solidus (% of summer)','FontSize',12);

title({'PCM Melting-Point Sensitivity Sweep (4 kg per panel)';...
       sprintf('Jabrayil, Azerbaijan | JJA TMY | state-dependent h_{pvcm} (%.0f/%.0f/%.0f W/m²·K) | m_{pcm}=%.0f kg', ...
               h_s, h_m, h_l, m_pcm)},...
      'FontSize',11,'FontWeight','normal');

legend({'Energy gain (%)', sprintf('Least-unfavourable: T_{melt}=%.0f°C', T_melt_range(best_idx)), ...
        'PCM at/above solidus (%)'},...
       'Location','northwest','FontSize',9);

xticks(T_melt_range);
grid on; box on;
text(T_melt_range(best_idx), E_gain_pct(best_idx)+ylim_left*0.12, ...
     sprintf('Least-unfavourable\n%.0f°C', T_melt_range(best_idx)), ...
     'HorizontalAlignment','center','FontSize',9,'Color','red','FontWeight','bold');

fprintf('Sensitivity figure generated.\n');

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