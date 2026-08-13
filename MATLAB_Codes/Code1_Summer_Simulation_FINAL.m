% =========================================================
%  CODE 1: PV-PCM SUMMER SIMULATION (FINAL CONFIGURATION)
%  ---------------------------------------------------------
%  Author: FARID ALIYEV
%  Date:   August, 2026
%  ---------------------------------------------------------
%  Dissertation: Fatty Acid PCM Cooling of Photovoltaic
%  Panels — A Simulation Study for Karabakh, Azerbaijan
%
%  Simulates one JJA summer (2,208 hourly records) of a PV
%  panel with and without a rear-mounted fatty acid PCM
%  layer, using the two-node lumped-capacitance model with
%  enthalpy-method phase change (Sections 3.6–3.7 of the
%  report). Figures use local time (UTC+4); PVGIS
%  timestamps are UTC.
%
%  CONFIGURATION:
%    m_pcm = 30 kg per panel (18.3 kg/m²)
%    Case A: MA-SA 64:36     (T_melt 44.13°C, L 182.4 J/g)
%    Case B: LA-SA 75.5:24.5 (T_melt 37.0°C,  L 182.7 J/g)
%
%  HOW TO RUN:
%    1. Place this script and the input files below in the
%       same working directory (MATLAB R2020a or later).
%    2. Run:  Code1_Summer_Simulation_FINAL
%       (runtime one to two minutes on a desktop PC)
%
%  INPUTS:  PVGIS_TMY_Summer_JJA.mat  (required)
%           sensitivity_results.csv and
%           mass_sensitivity_results.csv (outputs of the two
%           sensitivity scripts; needed only for Figures
%           4.5 and 4.6)
%  OUTPUTS: final_hourly_results.csv, final_summary_stats.csv
%           Figure_4_1_FINAL.png ... Figure_4_6_FINAL.png
%           (exported at 300 dpi)
% =========================================================

clc; clear; close all;

%% ── 1. LOAD DATA ─────────────────────────────────────────
DATA_PATH = 'PVGIS_TMY_Summer_JJA.mat';
fprintf('Loading PVGIS TMY summer data...\n');
d       = load(DATA_PATH);
G       = double(d.GHI(1,:));
T_amb   = double(d.T_amb(1,:));
v_wind  = double(d.wind_speed(1,:));
month   = double(d.month(1,:));
hour_utc= double(d.hour(1,:));
N       = length(G);
t_days  = (0:N-1) / 24;
hour_loc= mod(hour_utc + 4, 24);        % local time (UTC+4)
fprintf('  Loaded %d hourly records.\n\n', N);

%% ── 2. PV PANEL PARAMETERS (Table 3.3) ───────────────────
A_pv=1.638; m_pv=18.0; Cp_pv=900; alpha=0.90; eps_pv=0.91;
eta_ref=0.156; beta_T=0.0045; T_ref=25.0; sigma=5.67e-8;

%% ── 3. PCM PARAMETERS — FINAL CONFIGURATION ──────────────
m_pcm  = 30.0;        % PCM mass per panel [kg]
dT_melt= 4.0;
h_s=17; h_m=200; h_l=60;   % state-dependent conductance [W/m²·K]

% Case A: MA-SA 64:36 (Sarı & Kaygusuz, 2006)
TmA = 44.13;  LA_ = 182400;  CpsA = 2090;  CplA = 2200;
% Case B: LA-SA 75.5:24.5 (Sarı & Kaygusuz, 2002; Table 2.1)
TmB = 37.0;   LB_ = 182700;  CpsB = 2090;  CplB = 2200;

dt = 10;  n_sub = 3600/dt;

%% ── 4. BASELINE SIMULATION ───────────────────────────────
fprintf('Running BASELINE...\n');
[T_pv_B, P_B] = run_baseline(G,T_amb,v_wind,N,n_sub,dt, ...
    A_pv,m_pv,Cp_pv,alpha,eps_pv,eta_ref,beta_T,T_ref,sigma);

%% ── 5. PCM SIMULATIONS ───────────────────────────────────
fprintf('Running PCM case A (MA-SA, 30 kg)...\n');
[T_pv_A, T_pcm_Aa, P_A] = run_pcm(G,T_amb,v_wind,N,n_sub,dt, ...
    A_pv,m_pv,Cp_pv,alpha,eps_pv,eta_ref,beta_T,T_ref,sigma, ...
    TmA,dT_melt,LA_,CpsA,CplA,m_pcm,h_s,h_m,h_l);

fprintf('Running PCM case B (LA-SA, 30 kg)...\n');
[T_pv_Bb, T_pcm_Bb, P_Bb] = run_pcm(G,T_amb,v_wind,N,n_sub,dt, ...
    A_pv,m_pv,Cp_pv,alpha,eps_pv,eta_ref,beta_T,T_ref,sigma, ...
    TmB,dT_melt,LB_,CpsB,CplB,m_pcm,h_s,h_m,h_l);

%% ── 6. RESULTS ───────────────────────────────────────────
eta_B = max(eta_ref*(1-beta_T*(T_pv_B -T_ref)),0);
eta_A = max(eta_ref*(1-beta_T*(T_pv_A -T_ref)),0);
eta_C = max(eta_ref*(1-beta_T*(T_pv_Bb-T_ref)),0);

E_B = sum(P_B)/1000;  E_A = sum(P_A)/1000;  E_C = sum(P_Bb)/1000;
gA = E_A-E_B; gC = E_C-E_B;
idx_gen = G>10;

fprintf('\n══════════════════════════════════════════════════════\n');
fprintf('  FINAL RESULTS — Jabrayil JJA TMY | m_pcm = %.0f kg\n', m_pcm);
fprintf('══════════════════════════════════════════════════════\n');
fprintf('  Baseline energy:         %.4f kWh\n', E_B);
fprintf('  MA-SA (44.13°C) energy:  %.4f kWh  | gain %+.4f kWh (%+.4f%%)\n', E_A, gA, gA/E_B*100);
fprintf('  LA-SA (37.0°C)  energy:  %.4f kWh  | gain %+.4f kWh (%+.4f%%)\n', E_C, gC, gC/E_B*100);
fprintf('  Peak T:  baseline %.2f | MA-SA %.2f | LA-SA %.2f °C\n', max(T_pv_B), max(T_pv_A), max(T_pv_Bb));
fprintf('  Mean daytime T: baseline %.2f | MA-SA %.2f | LA-SA %.2f °C\n', ...
        mean(T_pv_B(idx_gen)), mean(T_pv_A(idx_gen)), mean(T_pv_Bb(idx_gen)));
fprintf('  Mean daytime dT: MA-SA %+.2f | LA-SA %+.2f °C\n', ...
        mean(T_pv_B(idx_gen)-T_pv_A(idx_gen)), mean(T_pv_B(idx_gen)-T_pv_Bb(idx_gen)));
fprintf('  PCM re-solidification (min T): MA-SA %.1f | LA-SA %.1f °C\n\n', min(T_pcm_Aa), min(T_pcm_Bb));

% Monthly
months_jja=[6 7 8]; month_names={'June','July','August'};
for m=1:3
    idx = month==months_jja(m);
    EmB(m)=sum(P_B(idx))/1000; EmA(m)=sum(P_A(idx))/1000; EmC(m)=sum(P_Bb(idx))/1000;
end

%% ── 7. EXPORT CSVs ───────────────────────────────────────
fid=fopen('final_hourly_results.csv','w');
fprintf(fid,'hour_index,month,hour_utc,hour_local,G_Wm2,T_amb_C,T_pv_baseline_C,T_pv_MASA_C,T_pv_LASA_C,T_pcm_MASA_C,T_pcm_LASA_C,P_baseline_W,P_MASA_W,P_LASA_W\n');
for i=1:N
    fprintf(fid,'%d,%d,%d,%d,%.2f,%.2f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n', ...
        i,month(i),hour_utc(i),hour_loc(i),G(i),T_amb(i),T_pv_B(i),T_pv_A(i),T_pv_Bb(i), ...
        T_pcm_Aa(i),T_pcm_Bb(i),P_B(i),P_A(i),P_Bb(i));
end
fclose(fid);

fid=fopen('final_summary_stats.csv','w');
fprintf(fid,'metric,June,July,August,JJA_total\n');
fprintf(fid,'Baseline_energy_kWh,%.4f,%.4f,%.4f,%.4f\n',EmB,E_B);
fprintf(fid,'MASA_energy_kWh,%.4f,%.4f,%.4f,%.4f\n',EmA,E_A);
fprintf(fid,'LASA_energy_kWh,%.4f,%.4f,%.4f,%.4f\n',EmC,E_C);
fprintf(fid,'MASA_gain_kWh,%.4f,%.4f,%.4f,%.4f\n',EmA-EmB,gA);
fprintf(fid,'LASA_gain_kWh,%.4f,%.4f,%.4f,%.4f\n',EmC-EmB,gC);
fprintf(fid,'MASA_gain_pct,%.4f,%.4f,%.4f,%.4f\n',(EmA-EmB)./EmB*100,gA/E_B*100);
fprintf(fid,'LASA_gain_pct,%.4f,%.4f,%.4f,%.4f\n',(EmC-EmB)./EmB*100,gC/E_B*100);
for m=1:3
    idx=month==months_jja(m);
    pkB(m)=max(T_pv_B(idx)); pkA(m)=max(T_pv_A(idx)); pkC(m)=max(T_pv_Bb(idx));
end
fprintf(fid,'Peak_T_baseline_C,%.2f,%.2f,%.2f,%.2f\n',pkB,max(T_pv_B));
fprintf(fid,'Peak_T_MASA_C,%.2f,%.2f,%.2f,%.2f\n',pkA,max(T_pv_A));
fprintf(fid,'Peak_T_LASA_C,%.2f,%.2f,%.2f,%.2f\n',pkC,max(T_pv_Bb));
fprintf(fid,'Min_T_pcm_MASA_C,%.2f,%.2f,%.2f,%.2f\n',...
    min(T_pcm_Aa(month==6)),min(T_pcm_Aa(month==7)),min(T_pcm_Aa(month==8)),min(T_pcm_Aa));
fprintf(fid,'Min_T_pcm_LASA_C,%.2f,%.2f,%.2f,%.2f\n',...
    min(T_pcm_Bb(month==6)),min(T_pcm_Bb(month==7)),min(T_pcm_Bb(month==8)),min(T_pcm_Bb));
fclose(fid);
fprintf('CSVs saved.\n');

%% ── 8. FIGURES (dissertation numbering; titles descriptive) ──
set(0,'DefaultAxesFontName','Times New Roman','DefaultAxesFontSize',11, ...
      'DefaultTextFontName','Times New Roman','DefaultLineLineWidth',1.5);
C_bl=[0.85 0.15 0.15]; C_ma=[0.12 0.47 0.71]; C_la=[0.20 0.63 0.17]; C_am=[0.6 0.6 0.6];

% Figure 4.1 — (a) daily maxima, full summer; (b) sustained 3-day hot episode
ndays   = N/24;
dayix   = repelem((1:ndays)',24);
dmax_B  = accumarray(dayix, T_pv_B(:),  [], @max);
dmax_A  = accumarray(dayix, T_pv_A(:),  [], @max);
dmax_am = accumarray(dayix, T_amb(:),   [], @max);
m3 = movmean(dmax_B,[2 0]);            % trailing 3-day mean of daily maxima
[~,we] = max(m3); ws = we-2;           % sustained-heat window (1-based days)
hsel = (dayix>=ws & dayix<=we);

f1=figure('Position',[50 50 900 640],'Color','w');
subplot(2,1,1); hold on;
plot(0:ndays-1,dmax_am,'-','Color',C_am,'LineWidth',1.2,'DisplayName','Ambient (daily max)');
plot(0:ndays-1,dmax_B ,'-','Color',C_bl,'LineWidth',1.7,'DisplayName','Baseline (daily max)');
plot(0:ndays-1,dmax_A ,'-','Color',C_ma,'LineWidth',1.7,'DisplayName','MA-SA PCM 30 kg (daily max)');
yline(TmA,':k','LineWidth',1.1,'DisplayName','Melting point (44.1°C)');
xlabel('Day of simulation (1 June = Day 0)'); ylabel('Daily max temperature (°C)');
title('(a) Daily maximum panel temperature, full summer','FontSize',10,'FontWeight','normal');
legend('Location','southwest','FontSize',8,'NumColumns',2);
grid on; box on; xlim([0 ndays-1]); set(gca,'GridAlpha',0.3);

subplot(2,1,2); hold on;
tt=(0:sum(hsel)-1)/24;
plot(tt,T_amb(hsel)   ,'-','Color',C_am,'LineWidth',1.1,'DisplayName','Ambient');
plot(tt,T_pv_B(hsel)  ,'-','Color',C_bl,'LineWidth',1.7,'DisplayName','Baseline');
plot(tt,T_pv_A(hsel)  ,'-','Color',C_ma,'LineWidth',1.7,'DisplayName','MA-SA PCM 30 kg');
plot(tt,T_pcm_Aa(hsel),'--','Color',C_la,'LineWidth',1.3,'DisplayName','PCM temperature');
yline(TmA,':k','LineWidth',1.1,'DisplayName','Melting point (44.1°C)');
xlabel(sprintf('Day within sustained hot episode (Days %d–%d)',ws-1,we-1));
ylabel('Temperature (°C)');
title('(b) Sustained 3-day hot episode, hourly detail','FontSize',10,'FontWeight','normal');
legend('Location','northwest','FontSize',8,'NumColumns',2);
grid on; box on; xlim([0 3]); set(gca,'GridAlpha',0.3);
exportgraphics(f1,'Figure_4_1_FINAL.png','Resolution',300);

% Figure 4.2 — efficiency over an average July day (local time)
eB_h=zeros(1,24); eA_h=zeros(1,24); g_h=zeros(1,24);
idx_jul2 = month==7;
for h=0:23
    ih = idx_jul2 & (hour_loc==h);
    eB_h(h+1)=mean(eta_B(ih))*100; eA_h(h+1)=mean(eta_A(ih))*100;
    g_h(h+1)=mean(G(ih));
end
f2=figure('Position',[50 500 900 420],'Color','w'); hold on;
yl2=[12.8 16.8];
night = g_h<5;
dn=diff([0 night 0]); ns=find(dn==1)-1; ne=find(dn==-1)-2;
for k=1:numel(ns)
    patch([ns(k)-0.5 ne(k)+0.5 ne(k)+0.5 ns(k)-0.5], ...
          [yl2(1) yl2(1) yl2(2) yl2(2)],[0.88 0.88 0.88], ...
          'EdgeColor','none','HandleVisibility','off');
end
patch(NaN(1,4),NaN(1,4),[0.88 0.88 0.88],'EdgeColor','none', ...
      'DisplayName','Night (no generation)');
plot(0:23,eB_h,'-','Color',C_bl,'LineWidth',1.9,'DisplayName','Baseline');
plot(0:23,eA_h,'-','Color',C_ma,'LineWidth',1.9,'DisplayName','MA-SA PCM 30 kg');
yline(15.6,':k','LineWidth',1.1,'DisplayName','Reference 15.6% (25°C)');
xlabel('Hour of day (local time, UTC+4)'); ylabel('Mean PV efficiency (%)');
title('PV Efficiency over an Average July Day (Jabrayil, TMY)','FontSize',10,'FontWeight','normal');
legend('Location','southwest','FontSize',8); grid on; box on;
xlim([0 23]); ylim(yl2); xticks(0:2:23); set(gca,'GridAlpha',0.3);
exportgraphics(f2,'Figure_4_2_FINAL.png','Resolution',300);

% Figure 4.3 — monthly energy, three cases
f3=figure('Position',[950 50 640 420],'Color','w');
bd=[EmB; EmA; EmC]';
b=bar(bd,'grouped'); b(1).FaceColor=C_bl; b(2).FaceColor=C_ma; b(3).FaceColor=C_la;
hold on;
for m=1:3
    text(m, max(bd(m,:))+1.2, sprintf('MA-SA %+.2f%%  LA-SA %+.2f%%', ...
        (EmA(m)-EmB(m))/EmB(m)*100, (EmC(m)-EmB(m))/EmB(m)*100), ...
        'FontSize',8,'HorizontalAlignment','center');
end
xlabel('Month'); ylabel('Energy Generated per Panel (kWh)');
title('Monthly Energy Output: Baseline vs. MA-SA and LA-SA PCM, 30 kg (Jabrayil, JJA TMY)','FontSize',10,'FontWeight','normal');
legend({'Baseline (no PCM)','MA-SA 64:36 (44.1°C)','LA-SA 75.5:24.5 (37.0°C)'},'Location','southoutside','Orientation','horizontal','FontSize',9);
set(gca,'XTickLabel',month_names,'GridAlpha',0.3); grid on; box on; ylim([0 60]);
exportgraphics(f3,'Figure_4_3_FINAL.png','Resolution',300);

% Figure 4.4 — average July day temperature profile, LOCAL time
f4=figure('Position',[950 500 640 420],'Color','w');
idx_jul = month==7;
for h=0:23
    ih = idx_jul & (hour_loc==h);
    havg_B(h+1)=mean(T_pv_B(ih)); havg_A(h+1)=mean(T_pv_A(ih));
    havg_C(h+1)=mean(T_pv_Bb(ih)); havg_am(h+1)=mean(T_amb(ih));
end
hold on;
plot(0:23,havg_am,'-','Color',C_am,'LineWidth',1.2,'DisplayName','Ambient');
plot(0:23,havg_B,'-','Color',C_bl,'LineWidth',2.0,'DisplayName','PV — Baseline');
plot(0:23,havg_A,'-','Color',C_ma,'LineWidth',2.0,'DisplayName','PV — MA-SA PCM (30 kg)');
plot(0:23,havg_C,'-','Color',C_la,'LineWidth',2.0,'DisplayName','PV — LA-SA PCM (30 kg)');
yline(TmA,':','Color',C_ma,'LineWidth',1.1,'DisplayName','MA-SA melting point (44.1°C)');
yline(TmB,':','Color',C_la,'LineWidth',1.1,'DisplayName','LA-SA melting point (37.0°C)');
xlabel('Hour of day (local time, UTC+4)'); ylabel('Mean Temperature (°C)');
title('Panel Temperature over an Average July Day (Jabrayil, TMY)','FontSize',10,'FontWeight','normal');
legend('Location','northoutside','NumColumns',3,'FontSize',8); grid on; box on; xlim([0 23]); xticks(0:2:23);
set(gca,'GridAlpha',0.3);
exportgraphics(f4,'Figure_4_4_FINAL.png','Resolution',300);

% Figure 4.5 — melting-point sweep (from sensitivity_results.csv)
try
    S = readmatrix('sensitivity_results.csv');
    f5=figure('Position',[100 100 980 440],'Color','w');
    yyaxis left
    bar(S(:,1),S(:,5),0.5,'FaceColor',C_ma,'FaceAlpha',0.85); hold on;
    [~,bi]=max(S(:,5));
    plot(S(bi,1),S(bi,5),'o','MarkerSize',9,'LineWidth',1.2, ...
        'MarkerEdgeColor',[0.06 0.35 0.15],'MarkerFaceColor',[0.13 0.55 0.24]);
    yline(0,'k-','LineWidth',1.0,'HandleVisibility','off');
    xlabel('PCM Melting Temperature (°C)'); ylabel('Summer Energy Gain (%)');
    yl=max(abs(S(:,5)))*1.4; ylim([-yl yl]);
    cand=[34.2 37.0 44.13]; clab={'LA-MA (34.2°C)','LA-SA (37°C)','MA-SA (44.1°C)'};
    for k=1:3
        xline(cand(k),':','Color',[0.3 0.3 0.3],'LineWidth',1.0,'HandleVisibility','off');
        text(cand(k),-yl*0.85,clab{k},'Rotation',90,'FontSize',8,'Color',[0.3 0.3 0.3]);
    end
    yyaxis right
    plot(S(:,1),S(:,7),'ko-','LineWidth',1.4,'MarkerFaceColor','k','MarkerSize',5);
    ylabel('Hours with PCM at/above Solidus (% of summer)');
    title('PCM Melting-Point Sensitivity at 4 kg (Jabrayil, JJA TMY)','FontSize',10,'FontWeight','normal');
    legend({'Energy gain (%)','Least-unfavourable (40°C)','PCM at/above solidus (%)'},'Location','northwest','FontSize',9);
    xticks(S(:,1)); grid on; box on;
    exportgraphics(f5,'Figure_4_5_FINAL.png','Resolution',300);
catch, fprintf('NOTE: sensitivity_results.csv not found — Figure 4.5 skipped\n'); end

% Figure 4.6 — mass sweep (from mass_sensitivity_results.csv)
try
    M = readmatrix('mass_sensitivity_results.csv');
    f6=figure('Position',[100 100 900 420],'Color','w'); hold on;
    i40 = abs(M(:,1)-40.0)<0.01;  iMA = abs(M(:,1)-44.13)<0.01;
    plot(M(i40,2),M(i40,7),'o-','Color',C_ma,'LineWidth',1.8,'MarkerFaceColor',C_ma,'DisplayName','T_{melt} = 40°C (sweep optimum)');
    plot(M(iMA,2),M(iMA,7),'s-','Color',C_bl,'LineWidth',1.8,'MarkerFaceColor',C_bl,'DisplayName','T_{melt} = 44.13°C (MA-SA)');
    yline(0,'k-','LineWidth',1.0,'HandleVisibility','off');
    xline(30,':','Color',[0.3 0.3 0.3],'LineWidth',1.1,'HandleVisibility','off');
    text(30.4,-0.1,'Final configuration (30 kg)','FontSize',8,'Color',[0.3 0.3 0.3],'Rotation',90);
    xlabel('PCM Mass per Panel (kg)'); ylabel('Summer Energy Gain (%)');
    title('PCM Mass Sensitivity (Jabrayil, JJA TMY)','FontSize',10,'FontWeight','normal');
    legend('Location','northwest','FontSize',10); grid on; box on; xticks(unique(M(:,2)));
    exportgraphics(f6,'Figure_4_6_FINAL.png','Resolution',300);
catch, fprintf('NOTE: mass_sensitivity_results.csv not found — Figure 4.6 skipped\n'); end

fprintf('All figures exported at 300 dpi.\nDONE - FINAL RUN COMPLETE\n');

%% ── LOCAL FUNCTIONS ──────────────────────────────────────
function [T_pv,P] = run_baseline(G,T_amb,v_wind,N,n_sub,dt,A_pv,m_pv,Cp_pv,alpha,eps_pv,eta_ref,beta_T,T_ref,sigma)
    T_pv=zeros(1,N); P=zeros(1,N); Tn=T_amb(1);
    for i=1:N
        Gi=G(i); Ta=T_amb(i); v=v_wind(i);
        hf=5.7+3.8*v; hb=0.5*(5.7+3.8*v);
        Tsky=0.0552*(Ta+273.15)^1.5;
        for s=1:n_sub
            TK=Tn+273.15;
            eta=max(eta_ref*(1-beta_T*(Tn-T_ref)),0);
            dT=(alpha*Gi*A_pv-eta*Gi*A_pv-hf*A_pv*(Tn-Ta) ...
               -eps_pv*sigma*A_pv*(TK^4-Tsky^4)-hb*A_pv*(Tn-Ta))/(m_pv*Cp_pv);
            Tn=Tn+dT*dt;
        end
        T_pv(i)=Tn;
        P(i)=max(eta_ref*(1-beta_T*(Tn-T_ref)),0)*Gi*A_pv;
    end
end

function [T_pv,T_pcm,P] = run_pcm(G,T_amb,v_wind,N,n_sub,dt,A_pv,m_pv,Cp_pv,alpha,eps_pv,eta_ref,beta_T,T_ref,sigma,Tm,dTm,L,Cps,Cpl,m_pcm,h_s,h_m,h_l)
    Tsol=Tm-dTm/2; Tliq=Tm+dTm/2;
    T_pv=zeros(1,N); T_pcm=zeros(1,N); P=zeros(1,N);
    Tn=T_amb(1); Tp=T_amb(1);
    H=pcm_T_to_H(Tp,Cps,Cpl,L,Tsol,Tliq);
    for i=1:N
        Gi=G(i); Ta=T_amb(i); v=v_wind(i);
        hf=5.7+3.8*v; hb=0.5*(5.7+3.8*v);
        Tsky=0.0552*(Ta+273.15)^1.5;
        for s=1:n_sub
            TK=Tn+273.15;
            eta=max(eta_ref*(1-beta_T*(Tn-T_ref)),0);
            if Tp<Tsol, he=h_s; elseif Tp<=Tliq, he=h_m; else, he=h_l; end
            Qp=he*A_pv*(Tn-Tp);
            dT=(alpha*Gi*A_pv-eta*Gi*A_pv-hf*A_pv*(Tn-Ta) ...
               -eps_pv*sigma*A_pv*(TK^4-Tsky^4)-Qp)/(m_pv*Cp_pv);
            Tn=Tn+dT*dt;
            H=H+(Qp-hb*A_pv*(Tp-Ta))/m_pcm*dt;
            Tp=pcm_H_to_T(H,Cps,Cpl,L,Tsol,Tliq);
        end
        T_pv(i)=Tn; T_pcm(i)=Tp;
        P(i)=max(eta_ref*(1-beta_T*(Tn-T_ref)),0)*Gi*A_pv;
    end
end

function H = pcm_T_to_H(T,Cps,Cpl,L,Tsol,Tliq)
    if T<Tsol, H=Cps*T;
    elseif T<=Tliq, H=Cps*Tsol+L*(T-Tsol)/(Tliq-Tsol);
    else, H=Cps*Tsol+L+Cpl*(T-Tliq); end
end
function T = pcm_H_to_T(H,Cps,Cpl,L,Tsol,Tliq)
    Hs=Cps*Tsol; Hl=Hs+L;
    if H<Hs, T=H/Cps;
    elseif H<=Hl, T=Tsol+(H-Hs)/L*(Tliq-Tsol);
    else, T=Tliq+(H-Hl)/Cpl; end
end
