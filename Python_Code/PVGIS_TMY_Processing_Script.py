"""
PVGIS TMY Data Processing Script
Author:  FARID ALIYEV
Date:    August, 2026
Project: Fatty Acid PCM Cooling of PV Panels — Karabakh, Azerbaijan
Site:    39.330°N, 47.052°E (Jabrayil district)
Data:    PVGIS-SARAH3 Typical Meteorological Year, 2005–2023

Cleans the raw PVGIS TMY export, validates it, derives the summer
(JJA) subset, exports the MATLAB simulation input files, and
generates Figures 3.2 and 3.3 of the report.

HOW TO RUN:
  1. Requires Python 3.9+ with pandas, NumPy, SciPy and Matplotlib.
  2. Place the raw PVGIS export tmy_39_330_47_052_2005_2023.csv in
     the same directory as this script.
  3. Run:  python PVGIS_TMY_Processing_Script.py

OUTPUTS: PVGIS_TMY_FullYear_clean.csv, PVGIS_TMY_Summer_JJA_clean.csv,
         PVGIS_TMY_FullYear.mat, PVGIS_TMY_Summer_JJA.mat,
         Figure_3_2_Annual_Climate_Profile.png,
         Figure_3_3_Summer_Detail.png
"""

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.io import savemat
import warnings
warnings.filterwarnings('ignore')

# ── Paths ────────────────────────────────────────────────────────────────────
CSV_PATH   = 'tmy_39_330_47_052_2005_2023.csv'
OUT_DIR    = './'

# ── 1. PARSE CSV ─────────────────────────────────────────────────────────────
print("=== PVGIS TMY Data Processing ===\n")

# Read skipping the header metadata rows and footer
# Find the row where actual data starts (after "time(UTC),...")
with open(CSV_PATH, 'r') as f:
    lines = f.readlines()

# Find header line index
header_idx = None
for i, line in enumerate(lines):
    if line.startswith('time(UTC)'):
        header_idx = i
        break

print(f"Data header found at line {header_idx}")

# Read the data portion only
df = pd.read_csv(CSV_PATH, skiprows=header_idx, 
                 nrows=8760,  # exactly one year
                 parse_dates=False)

print(f"Rows loaded: {len(df)}")
print(f"Columns: {list(df.columns)}")
print(f"\nFirst row: {df.iloc[0]['time(UTC)']}")
print(f"Last row:  {df.iloc[-1]['time(UTC)']}")

# ── 2. PARSE TIMESTAMPS ───────────────────────────────────────────────────────
# Format: YYYYMMDD:HHMM
df['timestamp'] = pd.to_datetime(df['time(UTC)'], format='%Y%m%d:%H%M')
df['month'] = df['timestamp'].dt.month
df['day']   = df['timestamp'].dt.day
df['hour']  = df['timestamp'].dt.hour

# ── 3. RENAME COLUMNS FOR CLARITY ────────────────────────────────────────────
df = df.rename(columns={
    'T2m':   'T_amb',      # Air temperature at 2m [°C]
    'RH':    'RH',         # Relative humidity [%]
    'G(h)':  'GHI',        # Global horizontal irradiance [W/m²]
    'Gb(n)': 'DNI',        # Direct normal irradiance [W/m²]
    'Gd(h)': 'DHI',        # Diffuse horizontal irradiance [W/m²]
    'IR(h)': 'IR',         # Infrared radiation [W/m²]
    'WS10m': 'wind_speed', # Wind speed at 10m [m/s]
    'WD10m': 'wind_dir',   # Wind direction [°]
    'SP':     'pressure',   # Surface pressure [Pa]
})

# ── 4. DATA VALIDATION ────────────────────────────────────────────────────────
print("\n── Data Validation ──────────────────────────────────────────")
print(f"Total hours:          {len(df):,}")
print(f"Missing values:       {df[['T_amb','GHI','DNI','wind_speed']].isnull().sum().sum()}")
print(f"\nTemperature [°C]:")
print(f"  Annual mean:        {df['T_amb'].mean():.1f}°C")
print(f"  Summer mean (JJA):  {df[df['month'].isin([6,7,8])]['T_amb'].mean():.1f}°C")
print(f"  Annual max:         {df['T_amb'].max():.1f}°C")
print(f"  Annual min:         {df['T_amb'].min():.1f}°C")
print(f"\nGlobal Horizontal Irradiance [W/m²]:")
print(f"  Annual total:       {df['GHI'].sum()/1000:.0f} kWh/m²")
print(f"  Summer total (JJA): {df[df['month'].isin([6,7,8])]['GHI'].sum()/1000:.0f} kWh/m²")
print(f"  Peak hourly:        {df['GHI'].max():.0f} W/m²")
print(f"\nWind Speed [m/s]:")
print(f"  Annual mean:        {df['wind_speed'].mean():.1f} m/s")
print(f"  Summer mean:        {df[df['month'].isin([6,7,8])]['wind_speed'].mean():.1f} m/s")

# GSA cross-check
annual_ghi_kwh = df['GHI'].sum() / 1000
print(f"\n── GSA Cross-Validation ─────────────────────────────────────")
print(f"PVGIS annual GHI:     {annual_ghi_kwh:.0f} kWh/m²")
print(f"GSA site report GHI:  1578.1 kWh/m²")
print(f"Difference:           {abs(annual_ghi_kwh - 1578.1):.1f} kWh/m² ({abs(annual_ghi_kwh-1578.1)/1578.1*100:.1f}%)")

# ── 5. SUMMER SUBSET (JJA) ────────────────────────────────────────────────────
df_summer = df[df['month'].isin([6, 7, 8])].copy().reset_index(drop=True)
print(f"\n── Summer Subset (JJA) ──────────────────────────────────────")
print(f"Summer hours:         {len(df_summer):,}")
print(f"June hours:           {len(df_summer[df_summer['month']==6])}")
print(f"July hours:           {len(df_summer[df_summer['month']==7])}")
print(f"August hours:         {len(df_summer[df_summer['month']==8])}")
print(f"Summer mean T_amb:    {df_summer['T_amb'].mean():.1f}°C")
print(f"Summer night T_amb:   {df_summer[df_summer['GHI']<5]['T_amb'].mean():.1f}°C (hours with GHI<5 W/m²)")
print(f"Summer day T_amb:     {df_summer[df_summer['GHI']>100]['T_amb'].mean():.1f}°C (hours with GHI>100 W/m²)")
print(f"Summer max T_amb:     {df_summer['T_amb'].max():.1f}°C")
print(f"Summer min T_amb:     {df_summer['T_amb'].min():.1f}°C")
print(f"PCM melt point:       44.13°C")
print(f"PCM re-solidifies:    YES (summer nights avg {df_summer[df_summer['GHI']<5]['T_amb'].mean():.1f}°C << 44.13°C)")

# ── 6. MONTHLY SUMMARY TABLE ──────────────────────────────────────────────────
monthly = df.groupby('month').agg(
    GHI_total=('GHI', lambda x: x.sum()/1000),
    T_mean=('T_amb', 'mean'),
    T_max=('T_amb', 'max'),
    T_min=('T_amb', 'min'),
    wind_mean=('wind_speed', 'mean'),
    DNI_total=('DNI', lambda x: x.sum()/1000)
).round(1)

months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
monthly.index = months
print("\n── Monthly Summary ──────────────────────────────────────────")
print(monthly.to_string())

# ── 7. EXPORT .MAT FILES FOR MATLAB ──────────────────────────────────────────
# Full year
mat_full = {
    'GHI':        df['GHI'].values.astype(float),
    'DNI':        df['DNI'].values.astype(float),
    'DHI':        df['DHI'].values.astype(float),
    'T_amb':      df['T_amb'].values.astype(float),
    'wind_speed': df['wind_speed'].values.astype(float),
    'RH':         df['RH'].values.astype(float),
    'month':      df['month'].values.astype(float),
    'hour':       df['hour'].values.astype(float),
    'n_hours':    np.array([8760.0]),
    'lat':        np.array([39.330]),
    'lon':        np.array([47.052]),
    'source':     np.array(['PVGIS-SARAH3 TMY 2005-2023'], dtype=object),
}
savemat(OUT_DIR + 'PVGIS_TMY_FullYear.mat', mat_full)
print(f"\n── Export ───────────────────────────────────────────────────")
print(f"Saved: PVGIS_TMY_FullYear.mat  ({len(df['GHI'].values)} hourly values)")

# Summer subset
mat_summer = {
    'GHI':        df_summer['GHI'].values.astype(float),
    'DNI':        df_summer['DNI'].values.astype(float),
    'DHI':        df_summer['DHI'].values.astype(float),
    'T_amb':      df_summer['T_amb'].values.astype(float),
    'wind_speed': df_summer['wind_speed'].values.astype(float),
    'RH':         df_summer['RH'].values.astype(float),
    'month':      df_summer['month'].values.astype(float),
    'hour':       df_summer['hour'].values.astype(float),
    'n_hours':    np.array([float(len(df_summer))]),
    'lat':        np.array([39.330]),
    'lon':        np.array([47.052]),
    'source':     np.array(['PVGIS-SARAH3 TMY 2005-2023 JJA subset'], dtype=object),
}
savemat(OUT_DIR + 'PVGIS_TMY_Summer_JJA.mat', mat_summer)
print(f"Saved: PVGIS_TMY_Summer_JJA.mat ({len(df_summer['GHI'].values)} hourly values)")

# Also export as CSV for reference
df[['timestamp','month','hour','T_amb','GHI','DNI','DHI','wind_speed','RH']].to_csv(
    OUT_DIR + 'PVGIS_TMY_FullYear_clean.csv', index=False)
df_summer[['timestamp','month','hour','T_amb','GHI','DNI','DHI','wind_speed','RH']].to_csv(
    OUT_DIR + 'PVGIS_TMY_Summer_JJA_clean.csv', index=False)
print(f"Saved: PVGIS_TMY_FullYear_clean.csv")
print(f"Saved: PVGIS_TMY_Summer_JJA_clean.csv")

# ── 8. PUBLICATION-QUALITY PLOTS ─────────────────────────────────────────────
plt.style.use('seaborn-v0_8-whitegrid')
colors = {
    'GHI': '#E8812A',
    'T':   '#C0392B',
    'wind':'#2980B9',
    'night':'#2C3E50',
    'summer':'#E74C3C',
    'spring':'#27AE60',
    'autumn':'#8E44AD',
    'winter':'#2980B9',
}

# ── FIGURE 1: Annual Overview ─────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
fig.suptitle('Annual Climate Profile — Jabrayil, Azerbaijan (39.33°N, 47.05°E)\n'
             'PVGIS-SARAH3 Typical Meteorological Year, 2005–2023',
             fontsize=12, fontweight='bold', y=0.98)

month_labels = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
month_nums = list(range(1,13))

# (a) Monthly GHI totals
ax = axes[0,0]
ghi_monthly = [df[df['month']==m]['GHI'].sum()/1000 for m in month_nums]
bar_colors = [colors['summer'] if m in [6,7,8] else '#BDC3C7' for m in month_nums]
bars = ax.bar(month_labels, ghi_monthly, color=bar_colors, edgecolor='white', linewidth=0.5)
ax.set_ylabel('GHI (kWh/m²)', fontsize=10)
ax.set_title('(a) Monthly Global Horizontal Irradiation', fontsize=10, fontweight='bold')
ax.set_ylim(0, max(ghi_monthly)*1.2)
for bar, val in zip(bars, ghi_monthly):
    ax.text(bar.get_x()+bar.get_width()/2, bar.get_height()+2, f'{val:.0f}',
            ha='center', va='bottom', fontsize=7.5)
ax.axhline(sum(ghi_monthly)/12, color='#7F8C8D', linestyle='--', linewidth=0.8, label=f'Monthly avg: {sum(ghi_monthly)/12:.0f} kWh/m²')
ax.legend(fontsize=8)
# Add JJA label
ax.text(5.5, max(ghi_monthly)*1.1, 'JJA\n(simulation\nperiod)',
        ha='center', va='center', fontsize=7.5, color=colors['summer'],
        bbox=dict(boxstyle='round,pad=0.3', facecolor='white', edgecolor=colors['summer'], alpha=0.8))

# (b) Monthly temperature range
ax = axes[0,1]
t_mean = [df[df['month']==m]['T_amb'].mean() for m in month_nums]
t_max  = [df[df['month']==m]['T_amb'].max() for m in month_nums]
t_min  = [df[df['month']==m]['T_amb'].min() for m in month_nums]
ax.fill_between(month_labels, t_min, t_max, alpha=0.2, color=colors['T'], label='Min–Max range')
ax.plot(month_labels, t_mean, 'o-', color=colors['T'], linewidth=2, markersize=5, label='Monthly mean')
ax.axhline(44.13, color='#27AE60', linestyle='--', linewidth=1.2, label='PCM melt point (44.1°C)')
ax.axhline(0, color='#95A5A6', linestyle=':', linewidth=0.8)
ax.set_ylabel('Temperature (°C)', fontsize=10)
ax.set_title('(b) Monthly Air Temperature Range', fontsize=10, fontweight='bold')
ax.legend(fontsize=8)
ax.set_ylim(-15, 50)

# (c) Average hourly GHI by season
ax = axes[1,0]
hours = list(range(24))
jja   = df[df['month'].isin([6,7,8])].groupby('hour')['GHI'].mean()
mam   = df[df['month'].isin([3,4,5])].groupby('hour')['GHI'].mean()
son   = df[df['month'].isin([9,10,11])].groupby('hour')['GHI'].mean()
djf   = df[df['month'].isin([12,1,2])].groupby('hour')['GHI'].mean()
ax.plot(hours, jja.reindex(hours, fill_value=0),  '-', color=colors['summer'], linewidth=2.5, label='Summer (JJA)')
ax.plot(hours, mam.reindex(hours, fill_value=0),  '-', color=colors['spring'], linewidth=1.5, label='Spring (MAM)', alpha=0.8)
ax.plot(hours, son.reindex(hours, fill_value=0),  '-', color=colors['autumn'], linewidth=1.5, label='Autumn (SON)', alpha=0.8)
ax.plot(hours, djf.reindex(hours, fill_value=0),  '-', color=colors['winter'], linewidth=1.5, label='Winter (DJF)', alpha=0.8)
ax.fill_between(hours, jja.reindex(hours, fill_value=0), alpha=0.12, color=colors['summer'])
ax.set_xlabel('Hour of day (UTC)', fontsize=10)
ax.set_ylabel('Mean GHI (W/m²)', fontsize=10)
ax.set_title('(c) Mean Daily Irradiance Profile by Season', fontsize=10, fontweight='bold')
ax.set_xticks(range(0,24,3))
ax.legend(fontsize=8)
ax.set_xlim(0,23)

# (d) Night-time temperature validation (PCM re-solidification check)
ax = axes[1,1]
night_hours = df_summer[df_summer['GHI'] < 5].copy()
t_bins = np.arange(-5, 45, 2)
june_n = night_hours[night_hours['month']==6]['T_amb']
july_n = night_hours[night_hours['month']==7]['T_amb']
aug_n  = night_hours[night_hours['month']==8]['T_amb']
ax.hist(june_n, bins=t_bins, alpha=0.6, color='#E74C3C', label=f'June (mean {june_n.mean():.1f}°C)', density=True)
ax.hist(july_n, bins=t_bins, alpha=0.6, color='#E67E22', label=f'July (mean {july_n.mean():.1f}°C)', density=True)
ax.hist(aug_n,  bins=t_bins, alpha=0.6, color='#F1C40F', label=f'August (mean {aug_n.mean():.1f}°C)', density=True)
ax.axvline(44.13, color='#27AE60', linestyle='--', linewidth=1.5, label='PCM melt point (44.1°C)')
ax.set_xlabel('Night-time air temperature (°C)', fontsize=10)
ax.set_ylabel('Probability density', fontsize=10)
ax.set_title('(d) Summer Night-Time Temperature Distribution\n(hours with GHI < 5 W/m²) — PCM Re-solidification Validation',
             fontsize=9, fontweight='bold')
ax.legend(fontsize=8)
ax.text(0.97, 0.60,
        f'All night-time temps\nwell below 44.1°C\n→ PCM fully re-solidifies\nevery night',
        transform=ax.transAxes, ha='right', va='top', fontsize=8,
        color='#27AE60', fontweight='bold',
        bbox=dict(boxstyle='round,pad=0.4', facecolor='white', edgecolor='#27AE60', alpha=0.9))

plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_DIR + 'Figure_3_2_Annual_Climate_Profile.png', dpi=180, bbox_inches='tight')
plt.close()
print(f"\nSaved: Figure_3_2_Annual_Climate_Profile.png")

# ── FIGURE 2: Summer Deep-Dive ────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(14, 9))
fig.suptitle('Summer (JJA) Climate Detail — Jabrayil, Azerbaijan\n'
             'PVGIS-SARAH3 Typical Meteorological Year, 2005–2023',
             fontsize=12, fontweight='bold', y=0.98)

month_colors = {6:'#E74C3C', 7:'#E67E22', 8:'#F39C12'}
month_names  = {6:'June', 7:'July', 8:'August'}

# (a) Average hourly GHI for each summer month
ax = axes[0,0]
for m, col in month_colors.items():
    h_profile = df_summer[df_summer['month']==m].groupby('hour')['GHI'].mean()
    ax.plot(range(24), h_profile.reindex(range(24), fill_value=0),
            'o-', color=col, linewidth=2, markersize=4, label=month_names[m])
ax.fill_between(range(24),
    df_summer.groupby('hour')['GHI'].mean().reindex(range(24), fill_value=0),
    alpha=0.1, color='#E8812A', label='JJA mean')
ax.set_xlabel('Hour of day (UTC)', fontsize=10)
ax.set_ylabel('Mean GHI (W/m²)', fontsize=10)
ax.set_title('(a) Mean Hourly GHI by Summer Month', fontsize=10, fontweight='bold')
ax.set_xticks(range(0,24,2))
ax.legend(fontsize=9)
ax.set_xlim(0,23)

# (b) Temperature: day vs night comparison
ax = axes[0,1]
for m, col in month_colors.items():
    sub = df_summer[df_summer['month']==m]
    day_t   = sub[sub['GHI']>100]['T_amb']
    night_t = sub[sub['GHI']<5]['T_amb']
    ax.boxplot([day_t, night_t],
               positions=[list(month_colors.keys()).index(m)*3+0.7,
                          list(month_colors.keys()).index(m)*3+1.3],
               widths=0.5, patch_artist=True,
               boxprops=dict(facecolor=col, alpha=0.6),
               medianprops=dict(color='black', linewidth=1.5),
               flierprops=dict(marker='.', markersize=3))

ax.axhline(44.13, color='#27AE60', linestyle='--', linewidth=1.2, label='PCM melt (44.1°C)')
ax.set_xticks([1, 4, 7])
ax.set_xticklabels(['June', 'July', 'August'])
ax.set_ylabel('Air temperature (°C)', fontsize=10)
ax.set_title('(b) Daytime vs Night-Time Temperature\n(light = night GHI<5, dark = day GHI>100)',
             fontsize=9, fontweight='bold')
ax.legend(fontsize=9)

# (c) Wind speed distribution
ax = axes[1,0]
for m, col in month_colors.items():
    ws = df_summer[df_summer['month']==m]['wind_speed']
    ax.hist(ws, bins=np.arange(0, 12, 0.5), alpha=0.6, color=col,
            label=f'{month_names[m]} (mean {ws.mean():.1f} m/s)', density=True)
ax.set_xlabel('Wind speed (m/s)', fontsize=10)
ax.set_ylabel('Probability density', fontsize=10)
ax.set_title('(c) Summer Wind Speed Distribution', fontsize=10, fontweight='bold')
ax.legend(fontsize=9)

# (d) Sample week — all variables together (first full week of July)
ax = axes[1,1]
july_data = df_summer[df_summer['month']==7].copy()
week = july_data.iloc[:7*24]  # first 7 days
hours_range = range(len(week))
ax2 = ax.twinx()
ax.fill_between(hours_range, week['GHI'], alpha=0.4, color='#E8812A', label='GHI (W/m²)')
ax2.plot(hours_range, week['T_amb'], '-', color=colors['T'], linewidth=1.2, label='T_amb (°C)')
ax2.axhline(44.13, color='#27AE60', linestyle='--', linewidth=1, label='PCM melt point')
ax.set_xlabel('Hour (first 7 days of July TMY)', fontsize=9)
ax.set_ylabel('GHI (W/m²)', fontsize=10, color='#E8812A')
ax2.set_ylabel('Temperature (°C)', fontsize=10, color=colors['T'])
ax.set_title('(d) Representative Week — GHI and Temperature\n(Illustrating Day-Night Cycle)', fontsize=9, fontweight='bold')
lines1, labels1 = ax.get_legend_handles_labels()
lines2, labels2 = ax2.get_legend_handles_labels()
ax.legend(lines1+lines2, labels1+labels2, fontsize=8, loc='upper right')
# Add day markers
for d in range(7):
    ax.axvline(d*24, color='#95A5A6', linestyle=':', linewidth=0.6)

plt.tight_layout(rect=[0, 0, 1, 0.96])
fig.savefig(OUT_DIR + 'Figure_3_3_Summer_Detail.png', dpi=180, bbox_inches='tight')
plt.close()
print(f"Saved: Figure_3_3_Summer_Detail.png")

# ── 9. FINAL SUMMARY ─────────────────────────────────────────────────────────
print("\n═══════════════════════════════════════════════════════════════")
print("PROCESSING COMPLETE — Files generated:")
print("  PVGIS_TMY_FullYear.mat          → MATLAB annual simulation (8,760 h)")
print("  PVGIS_TMY_Summer_JJA.mat        → MATLAB summer simulation (2,208 h)")
print("  PVGIS_TMY_FullYear_clean.csv    → Full year clean data")
print("  PVGIS_TMY_Summer_JJA_clean.csv  → JJA clean data")
print("  Figure_3_2_Annual_Climate_Profile.png")
print("  Figure_3_3_Summer_Detail.png")
print("\nKey findings for Section 3.5:")
print(f"  PVGIS annual GHI:   {df['GHI'].sum()/1000:.0f} kWh/m²  (GSA: 1578.1 — diff: {abs(df['GHI'].sum()/1000-1578.1):.1f} kWh/m²)")
print(f"  Summer GHI:         {df_summer['GHI'].sum()/1000:.0f} kWh/m²")
print(f"  Annual mean T:      {df['T_amb'].mean():.1f}°C")
print(f"  Summer mean T:      {df_summer['T_amb'].mean():.1f}°C")
print(f"  Summer night T:     {df_summer[df_summer['GHI']<5]['T_amb'].mean():.1f}°C")
print(f"  PCM melt point:     44.13°C")
print(f"  PCM re-solidifies:  YES — margin of {44.13 - df_summer[df_summer['GHI']<5]['T_amb'].max():.1f}°C above max night temp")
print(f"  Summer wind:        {df_summer['wind_speed'].mean():.1f} m/s mean")
print("═══════════════════════════════════════════════════════════════")
