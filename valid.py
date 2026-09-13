import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d

# =============================================================================
# 1. Leitura do arquivo experimental + CFAST
# =============================================================================
exp_file = "1Room1.csv"
df_exp = pd.read_csv(exp_file, skiprows=2, header=None)
print(f"Experimental: {len(df_exp)} linhas, {df_exp.shape[1]} colunas")

# ---------- Índices das colunas experimentais (verifique se estão corretos) ----------
IDX_TIME = 86
IDX_HRR = 8
IDX_INTRFACE = 10
IDX_TEMPUP = 83
IDX_TEMPLOW = 72
IDX_O2 = 21
IDX_CO2 = 2

# ---------- Extração dos dados experimentais ----------
t_exp = pd.to_numeric(df_exp[IDX_TIME], errors='coerce').values
hrr_exp = pd.to_numeric(df_exp[IDX_HRR], errors='coerce').values
z_int_exp = pd.to_numeric(df_exp[IDX_INTRFACE], errors='coerce').values
t_up_exp_raw = pd.to_numeric(df_exp[IDX_TEMPUP], errors='coerce').values
t_low_exp = pd.to_numeric(df_exp[IDX_TEMPLOW], errors='coerce').values
o2_exp_raw = pd.to_numeric(df_exp[IDX_O2], errors='coerce').values
co2_exp_raw = pd.to_numeric(df_exp[IDX_CO2], errors='coerce').values

# ---------- Extração dos dados CFAST (últimas 6 colunas) ----------
cfast_cols = df_exp.iloc[:, -6:]
t_cfast_min = pd.to_numeric(cfast_cols.iloc[:, 0], errors='coerce').values
t_up_cfast = pd.to_numeric(cfast_cols.iloc[:, 1], errors='coerce').values
z_cfast = pd.to_numeric(cfast_cols.iloc[:, 2], errors='coerce').values
o2_cfast_perc = pd.to_numeric(cfast_cols.iloc[:, 4], errors='coerce').values
co2_cfast_perc = pd.to_numeric(cfast_cols.iloc[:, 5], errors='coerce').values

t_cfast = t_cfast_min * 60.0
o2_cfast = o2_cfast_perc / 100.0
co2_cfast = co2_cfast_perc / 100.0

# Converter experimentais para fração e Celsius se necessário
if np.nanmax(o2_exp_raw) > 1.0:
    print("O2 experimental em porcentagem. Dividindo por 100.")
    o2_exp = o2_exp_raw / 100.0
else:
    o2_exp = o2_exp_raw

if np.nanmax(co2_exp_raw) > 1.0:
    print("CO2 experimental em porcentagem. Dividindo por 100.")
    co2_exp = co2_exp_raw / 100.0
else:
    co2_exp = co2_exp_raw

if np.nanmax(t_up_exp_raw) > 100:
    print("TEMPUP experimental em Kelvin. Convertendo para Celsius.")
    t_up_exp = t_up_exp_raw - 273.15
else:
    t_up_exp = t_up_exp_raw

# =============================================================================
# 2. Leitura das simulações (Heskestad, Zukoski, McCaffrey)
# =============================================================================
def load_simulacao(nome):
    # Lê arquivo, ignorando linhas que começam com '#'
    df = pd.read_csv(nome, comment='#')
    df.columns = df.columns.str.strip()
    # As primeiras colunas são: Tempo(s), Z_int(m), T_u(C), T_w(C), Y_O2, Y_CO2, ...
    # Pegamos as 6 primeiras
    t = df.iloc[:, 0].values
    z = df.iloc[:, 1].values
    T_up = df.iloc[:, 2].values
    T_low = df.iloc[:, 3].values
    o2 = df.iloc[:, 4].values
    co2 = df.iloc[:, 5].values
    return t, z, T_up, T_low, o2, co2

t_hes, z_hes, Tup_hes, Tlow_hes, o2_hes, co2_hes = load_simulacao("resultados_teste1_Babrauskas_heskestad.csv")
t_zuk, z_zuk, Tup_zuk, Tlow_zuk, o2_zuk, co2_zuk = load_simulacao("resultados_teste1_Babrauskas_zukoski.csv")
t_mcc, z_mcc, Tup_mcc, Tlow_mcc, o2_mcc, co2_mcc = load_simulacao("resultados_teste1_Babrauskas_mccaffrey.csv")

print(f"Heskestad: {len(t_hes)} pts, Zukoski: {len(t_zuk)} pts, McCaffrey: {len(t_mcc)} pts")

# =============================================================================
# 3. Interpolação para grid comum (0 a 2070 s)
# =============================================================================
t_grid = np.linspace(0, 2070, 2071)

def interp_clean(x, y, t_grid):
    mask = ~np.isnan(y) & ~np.isnan(x)
    if mask.sum() == 0:
        return np.full_like(t_grid, np.nan)
    f = interp1d(x[mask], y[mask], bounds_error=False, fill_value=np.nan)
    return f(t_grid)

# Experimentais
Tup_exp_i = interp_clean(t_exp, t_up_exp, t_grid)
Z_exp_i = interp_clean(t_exp, z_int_exp, t_grid)
HRR_exp_i = interp_clean(t_exp, hrr_exp, t_grid)
O2_exp_i = interp_clean(t_exp, o2_exp, t_grid)
CO2_exp_i = interp_clean(t_exp, co2_exp, t_grid)

# CFAST
Tup_cfast_i = interp_clean(t_cfast, t_up_cfast, t_grid)
Z_cfast_i = interp_clean(t_cfast, z_cfast, t_grid)
O2_cfast_i = interp_clean(t_cfast, o2_cfast, t_grid)
CO2_cfast_i = interp_clean(t_cfast, co2_cfast, t_grid)

# Simulações
Tup_hes_i = interp_clean(t_hes, Tup_hes, t_grid)
Z_hes_i = interp_clean(t_hes, z_hes, t_grid)
O2_hes_i = interp_clean(t_hes, o2_hes, t_grid)
CO2_hes_i = interp_clean(t_hes, co2_hes, t_grid)

Tup_zuk_i = interp_clean(t_zuk, Tup_zuk, t_grid)
Z_zuk_i = interp_clean(t_zuk, z_zuk, t_grid)
O2_zuk_i = interp_clean(t_zuk, o2_zuk, t_grid)
CO2_zuk_i = interp_clean(t_zuk, co2_zuk, t_grid)

Tup_mcc_i = interp_clean(t_mcc, Tup_mcc, t_grid)
Z_mcc_i = interp_clean(t_mcc, z_mcc, t_grid)
O2_mcc_i = interp_clean(t_mcc, o2_mcc, t_grid)
CO2_mcc_i = interp_clean(t_mcc, co2_mcc, t_grid)

# =============================================================================
# 4. Métricas de erro
# =============================================================================
def metrics(y_sim, y_ref, nome):
    mask = ~np.isnan(y_sim) & ~np.isnan(y_ref)
    if mask.sum() == 0:
        print(f"{nome}: sem dados")
        return
    y_s = y_sim[mask]
    y_r = y_ref[mask]
    rmse = np.sqrt(np.mean((y_s - y_r)**2))
    mae = np.mean(np.abs(y_s - y_r))
    mape = np.mean(np.abs((y_s - y_r) / (y_r + 1e-6))) * 100
    max_err = np.max(np.abs(y_s - y_r))
    print(f"{nome}: RMSE={rmse:.2f}, MAE={mae:.2f}, MAPE={mape:.2f}%, Max={max_err:.2f}")

print("\n=== MÉTRICAS  ===\n")
metrics(Tup_hes_i, Tup_exp_i, "Heskestad - T_up")
metrics(Tup_zuk_i, Tup_exp_i, "Zukoski - T_up")
metrics(Tup_mcc_i, Tup_exp_i, "McCaffrey - T_up")
metrics(Z_hes_i, Z_exp_i, "Heskestad - Z_int")
metrics(Z_zuk_i, Z_exp_i, "Zukoski - Z_int")
metrics(Z_mcc_i, Z_exp_i, "McCaffrey - Z_int")
metrics(O2_hes_i, O2_exp_i, "Heskestad - O2")
metrics(O2_zuk_i, O2_exp_i, "Zukoski - O2")
metrics(O2_mcc_i, O2_exp_i, "McCaffrey - O2")
metrics(CO2_hes_i, CO2_exp_i, "Heskestad - CO2")
metrics(CO2_zuk_i, CO2_exp_i, "Zukoski - CO2")
metrics(CO2_mcc_i, CO2_exp_i, "McCaffrey - CO2")

# =============================================================================
# 5. Gráficos
# =============================================================================
plt.style.use('seaborn-v0_8-darkgrid')

# 5.1 HRR (igual)
plt.figure(figsize=(10,5))
plt.plot(t_grid, HRR_exp_i, 'k-', label='Experimental (Babrauskas)')
plt.xlabel('Tempo (s)'); plt.ylabel('HRR (kW)')
plt.title('HRR – Teste 2')
plt.legend(); plt.grid(); plt.tight_layout()
plt.savefig('validacao_hrr.png', dpi=300); plt.close()

# 5.2 Temperatura – adiciona McCaffrey em verde
plt.figure(figsize=(10,5))
plt.plot(t_grid, Tup_cfast_i, 'ko', markersize=2, label='CFAST')
plt.plot(t_grid, Tup_exp_i, 'm.', markersize=1, label='Experimental')
plt.plot(t_hes, Tup_hes, 'r-', lw=1.5, label='Heskestad')
plt.plot(t_zuk, Tup_zuk, 'b-', lw=1.5, label='Zukoski')
plt.plot(t_mcc, Tup_mcc, 'g-', lw=1.5, label='McCaffrey')
plt.xlabel('Tempo (s)'); plt.ylabel('Temperatura (°C)')
plt.title('Temperatura da Camada Superior')
plt.legend(); plt.grid(); plt.tight_layout()
plt.savefig('validacao_temperatura.png', dpi=300); plt.close()

# 5.3 Interface – McCaffrey em verde
plt.figure(figsize=(10,5))
plt.plot(t_grid, Z_cfast_i, 'ko', markersize=2, label='CFAST')
plt.plot(t_grid, Z_exp_i, 'm.', markersize=1, label='Experimental')
plt.plot(t_hes, z_hes, 'r-', lw=1.5, label='Heskestad')
plt.plot(t_zuk, z_zuk, 'b-', lw=1.5, label='Zukoski')
plt.plot(t_mcc, z_mcc, 'g-', lw=1.5, label='McCaffrey')
plt.xlabel('Tempo (s)'); plt.ylabel('Altura da interface (m)')
plt.title('Altura da Interface')
plt.legend(); plt.grid(); plt.tight_layout()
plt.savefig('validacao_interface.png', dpi=300); plt.close()

# 5.4 O2 – McCaffrey em verde
plt.figure(figsize=(10,5))
plt.plot(t_grid, O2_cfast_i*100, 'ko', markersize=2, label='CFAST')
plt.plot(t_grid, O2_exp_i*100, 'm.', markersize=1, label='Experimental')
plt.plot(t_hes, o2_hes*100, 'r-', lw=1.5, label='Heskestad')
plt.plot(t_zuk, o2_zuk*100, 'b-', lw=1.5, label='Zukoski')
plt.plot(t_mcc, o2_mcc*100, 'g-', lw=1.5, label='McCaffrey')
plt.xlabel('Tempo (s)'); plt.ylabel('O₂ (%)')
plt.title('Fração mássica de Oxigênio')
plt.legend(); plt.grid(); plt.tight_layout()
plt.savefig('validacao_o2.png', dpi=300); plt.close()

# 5.5 CO2 – McCaffrey em verde
plt.figure(figsize=(10,5))
plt.plot(t_grid, CO2_cfast_i*100, 'ko', markersize=2, label='CFAST')
plt.plot(t_grid, CO2_exp_i*100, 'm.', markersize=1, label='Experimental')
plt.plot(t_hes, co2_hes*100, 'r-', lw=1.5, label='Heskestad')
plt.plot(t_zuk, co2_zuk*100, 'b-', lw=1.5, label='Zukoski')
plt.plot(t_mcc, co2_mcc*100, 'g-', lw=1.5, label='McCaffrey')
plt.xlabel('Tempo (s)'); plt.ylabel('CO₂ (%)')
plt.title('Fração mássica de Dióxido de Carbono')
plt.legend(); plt.grid(); plt.tight_layout()
plt.savefig('validacao_co2.png', dpi=300); plt.close()

# 5.6 Painel completo com todas as plumas
fig, axs = plt.subplots(5, 1, figsize=(10, 15), sharex=True)

# HRR
axs[0].plot(t_grid, HRR_exp_i, 'k-', label='Experimental')
axs[0].set_ylabel('HRR (kW)')
axs[0].legend()
axs[0].grid()

# Temperatura 
axs[1].plot(t_grid, Tup_cfast_i, 'ko', markersize=1, label='CFAST')
axs[1].plot(t_grid, Tup_exp_i, 'm.', markersize=1, label='Experimental')
axs[1].plot(t_hes, Tup_hes, 'r-', lw=1, label='Heskestad')
axs[1].plot(t_zuk, Tup_zuk, 'b-', lw=1, label='Zukoski')
axs[1].plot(t_mcc, Tup_mcc, 'g-', lw=1, label='McCaffrey')
axs[1].set_ylabel('T (°C)')
axs[1].legend()
axs[1].grid()

# Interface 
axs[2].plot(t_grid, Z_cfast_i, 'ko', markersize=1, label='CFAST')
axs[2].plot(t_grid, Z_exp_i, 'm.', markersize=1, label='Experimental')
axs[2].plot(t_hes, z_hes, 'r-', lw=1, label='Heskestad')
axs[2].plot(t_zuk, z_zuk, 'b-', lw=1, label='Zukoski')
axs[2].plot(t_mcc, z_mcc, 'g-', lw=1, label='McCaffrey')
axs[2].set_ylabel('Z_int (m)')
axs[2].legend()
axs[2].grid()

# O2 
axs[3].plot(t_grid, O2_cfast_i*100, 'ko', markersize=1, label='CFAST')
axs[3].plot(t_grid, O2_exp_i*100, 'm.', markersize=1, label='Experimental')
axs[3].plot(t_hes, o2_hes*100, 'r-', lw=1, label='Heskestad')
axs[3].plot(t_zuk, o2_zuk*100, 'b-', lw=1, label='Zukoski')
axs[3].plot(t_mcc, o2_mcc*100, 'g-', lw=1, label='McCaffrey')
axs[3].set_ylabel('O₂ (%)')
axs[3].legend()
axs[3].grid()

# CO2 
axs[4].plot(t_grid, CO2_cfast_i*100, 'ko', markersize=1, label='CFAST')
axs[4].plot(t_grid, CO2_exp_i*100, 'm.', markersize=1, label='Experimental')
axs[4].plot(t_hes, co2_hes*100, 'r-', lw=1, label='Heskestad')
axs[4].plot(t_zuk, co2_zuk*100, 'b-', lw=1, label='Zukoski')
axs[4].plot(t_mcc, co2_mcc*100, 'g-', lw=1, label='McCaffrey')
axs[4].set_xlabel('Tempo (s)')
axs[4].set_ylabel('CO₂ (%)')
axs[4].legend()
axs[4].grid()

fig.suptitle('Validação do Modelo – Teste 2 de Babrauskas (1984)', fontweight='bold')
fig.tight_layout()
plt.savefig('validacao_todas.png', dpi=300)
plt.close()

print("\nTodos os gráficos salvos como PNG.")
print("Validação concluída!")