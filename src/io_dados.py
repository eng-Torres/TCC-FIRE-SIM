"""
Leitura de arquivos de dados para validação.

Formatos suportados:
- simulador              : CSV do próprio simulador (resultados_*.csv)
- babrauskas_1984        : 1Room1.csv do NIST (experimental + CFAST)
- cfast_simple_compartments : Simple_compartments.csv (CFAST)
- cfast_simple_zone      : Simple_zone.csv (CFAST)

Cada leitor retorna um dict padronizado:
    {
        't':       ndarray [s],
        'T_up_C':  ndarray [°C]   ou None,
        'Z_int_m': ndarray [m]    ou None,
        'O2_vol':  ndarray [0-1]  ou None,
        'CO2_vol': ndarray [0-1]  ou None,
        'HRR_kW':  ndarray [kW]   ou None,
        'fonte':   str,
    }

Espécies são SEMPRE convertidas para fração volumétrica (0-1).
"""
from pathlib import Path

import numpy as np
import pandas as pd

# Massas molares [g/mol] — para conversão massa ↔ volume
M_AR  = 28.97
M_O2  = 32.00
M_CO2 = 44.01
M_CO  = 28.01


def _mass_para_vol(Y, M):
    """Fração mássica → fração volumétrica. X_i = Y_i · (M_ar / M_i)."""
    return Y * (M_AR / M)


# ---------------------------------------------------------------------------
# Leitor: CSV do próprio simulador
# ---------------------------------------------------------------------------
def ler_simulador(caminho, **_):
    """
    Lê CSV gerado por outputs.py.
    Colunas: Tempo(s), Z_int(m), T_u(C), T_w(C), Y_O2, Y_CO2, Y_CO, m_dot_in, m_dot_out.
    Converte Y_O2, Y_CO2 (mássica) para fração volumétrica.
    """
    df = pd.read_csv(caminho, comment='#')
    df.columns = [c.strip() for c in df.columns]

    t       = df.iloc[:, 0].to_numpy(float)
    Z_int   = df.iloc[:, 1].to_numpy(float)
    T_up_C  = df.iloc[:, 2].to_numpy(float)
    Y_O2    = df.iloc[:, 4].to_numpy(float)
    Y_CO2   = df.iloc[:, 5].to_numpy(float)

    return {
        't':       t,
        'T_up_C':  T_up_C,
        'Z_int_m': Z_int,
        'O2_vol':  _mass_para_vol(np.clip(Y_O2, 0.0, 1.0), M_O2),
        'CO2_vol': _mass_para_vol(np.clip(Y_CO2, 0.0, 1.0), M_CO2),
        'HRR_kW':  None,
        'fonte':   'simulador',
    }


# ---------------------------------------------------------------------------
# Leitor: Babrauskas 1984 (1Room1.csv do NIST)
# ---------------------------------------------------------------------------
# Índices das colunas experimentais (fixos no arquivo do NIST)
_BAB_IDX_TIME     = 86
_BAB_IDX_HRR      = 8
_BAB_IDX_INTERFACE = 10
_BAB_IDX_TEMPUP   = 83
_BAB_IDX_TEMPLOW  = 72
_BAB_IDX_O2       = 21
_BAB_IDX_CO2      = 2

# Últimas 6 colunas: NORMAL, Comp1(Time min), Comp1(?), Comp1(?), Comp1(O2 %), Comp1(CO2 %)
_BAB_CFAST_NCOLS = 6


def ler_babrauskas_1984(caminho, extrair_cfast=True, **_):
    """
    Lê o arquivo 1Room1.csv do NIST.

    Retorna o experimental como série principal. Se extrair_cfast=True,
    retorna também um dict 'cfast' no resultado.
    """
    df = pd.read_csv(caminho, skiprows=2, header=None)

    def _num(col):
        return pd.to_numeric(df[col], errors='coerce').to_numpy(float)

    t       = _num(_BAB_IDX_TIME)
    HRR     = _num(_BAB_IDX_HRR)/1000
    Z_int   = _num(_BAB_IDX_INTERFACE)
    T_up    = _num(_BAB_IDX_TEMPUP)
    O2_raw  = _num(_BAB_IDX_O2)
    CO2_raw = _num(_BAB_IDX_CO2)

    # Detecção de unidade (K vs °C) e (fração vs %)
    if np.nanmax(T_up) > 200:
        T_up = T_up - 273.15
    if np.nanmax(O2_raw) > 1.0:
        O2 = O2_raw / 100.0
    else:
        O2 = O2_raw
    if np.nanmax(CO2_raw) > 1.0:
        CO2 = CO2_raw / 100.0
    else:
        CO2 = CO2_raw

    resultado = {
        't':       t,
        'T_up_C':  T_up,
        'Z_int_m': Z_int,
        'O2_vol':  O2,
        'CO2_vol': CO2,
        'HRR_kW':  HRR,
        'fonte':   'babrauskas_1984 (experimental)',
    }

    if extrair_cfast and df.shape[1] >= _BAB_CFAST_NCOLS:
        cfast = df.iloc[:, -_BAB_CFAST_NCOLS:]

        t_min  = pd.to_numeric(cfast.iloc[:, 0], errors='coerce').to_numpy(float)
        T_cf   = pd.to_numeric(cfast.iloc[:, 1], errors='coerce').to_numpy(float)
        Z_cf   = pd.to_numeric(cfast.iloc[:, 2], errors='coerce').to_numpy(float)
        # índice 3 é Pressure — ignorado
        O2_cf  = pd.to_numeric(cfast.iloc[:, 4], errors='coerce').to_numpy(float)
        CO2_cf = pd.to_numeric(cfast.iloc[:, 5], errors='coerce').to_numpy(float)

        t_cf = t_min * 60.0     # min → s

        if np.nanmax(O2_cf)  > 1.0:
            O2_cf = O2_cf / 100.0
        if np.nanmax(CO2_cf) > 1.0:
            CO2_cf = CO2_cf / 100.0

        resultado['cfast'] = {
            't':       t_cf,
            'T_up_C':  T_cf,
            'Z_int_m': Z_cf,
            'O2_vol':  O2_cf,
            'CO2_vol': CO2_cf,
            'HRR_kW':  None,
            'fonte':   'babrauskas_1984 (CFAST)',
        }
    return resultado


# ---------------------------------------------------------------------------
# Leitor: CFAST Simple_compartments.csv
# ---------------------------------------------------------------------------
def ler_cfast_simple_compartments(caminho, **_):
    """
    Simple_compartments.csv:
    - Separador ';'
    - 4 linhas de cabeçalho: variáveis, descrição, compartimento, unidades
    - Colunas relevantes: Time, ULT_1, LLT_1, HGT_1, ULO2_1, ULCO2_1
    """
    df = pd.read_csv(caminho, sep=';', skiprows=4)
    df.columns = [c.strip() for c in df.columns]

    def _col(nome):
        if nome not in df.columns:
            return None
        return pd.to_numeric(df[nome], errors='coerce').to_numpy(float)

    t     = _col('Time')
    T_up  = _col('ULT_1')
    Z_int = _col('HGT_1')
    O2    = _col('ULO2_1')
    CO2   = _col('ULCO2_1')

    # CFAST reporta O2/CO2 em mol % (0–100) → converter para 0–1
    if O2  is not None and np.nanmax(O2)  > 1.0:
        O2 = O2 / 100.0
    if CO2 is not None and np.nanmax(CO2) > 1.0:
        CO2 = CO2 / 100.0

    return {
        't':       t,
        'T_up_C':  T_up,
        'Z_int_m': Z_int,
        'O2_vol':  O2,
        'CO2_vol': CO2,
        'HRR_kW':  None,
        'fonte':   'cfast_simple_compartments',
    }


# ---------------------------------------------------------------------------
# Leitor: CFAST Simple_zone.csv
# ---------------------------------------------------------------------------
def ler_cfast_simple_zone(caminho, **_):
    """
    Simple_zone.csv:
    - Separador ','
    - 2 linhas de cabeçalho: unidades, variáveis
    - Colunas relevantes: Time, ULT_1, LLT_1, HGT_1
    (Não há O2/CO2 neste arquivo.)
    """
    df = pd.read_csv(caminho, skiprows=2, sep=',', engine='python')
    df.columns = [c.strip() for c in df.columns]

    def _col(nome):
        if nome not in df.columns:
            return None
        return pd.to_numeric(df[nome], errors='coerce').to_numpy(float)

    return {
        't':       _col('Time'),
        'T_up_C':  _col('ULT_1'),
        'Z_int_m': _col('HGT_1'),
        'O2_vol':  None,
        'CO2_vol': None,
        'HRR_kW':  _col('HRR_1'),
        'fonte':   'cfast_simple_zone',
    }


# ---------------------------------------------------------------------------
# Registro de formatos + autodetecção + API pública
# ---------------------------------------------------------------------------
FORMATOS = {
    'simulador':                 ler_simulador,
    'babrauskas_1984':           ler_babrauskas_1984,
    'cfast_simple_compartments': ler_cfast_simple_compartments,
    'cfast_simple_zone':         ler_cfast_simple_zone,
}


def detectar_formato(caminho):
    """Detecta o formato do CSV pelas primeiras linhas."""
    caminho = Path(caminho)
    with open(caminho, encoding='utf-8', errors='ignore') as f:
        linhas = [f.readline() for _ in range(6)]

    l0 = linhas[0].lstrip()

    if l0.startswith('#'):
        return 'simulador'
    if l0.startswith('C1,C2,C3'):
        return 'babrauskas_1984'
    if l0.startswith('Time;ULT_1') or ';ULT_1;' in l0:
        return 'cfast_simple_compartments'
    if l0.startswith('s,C,C,m'):
        return 'cfast_simple_zone'
    if 'Time,ULT_1,LLT_1' in linhas[1] if len(linhas) > 1 else False:
        return 'cfast_simple_zone'

    return None


def ler_serie(caminho, formato=None, **kwargs):
    """
    Lê um arquivo de dados e retorna o dict padronizado.

    Se formato=None, tenta autodetectar pelo conteúdo.
    """
    caminho = Path(caminho)
    if not caminho.is_file():
        raise FileNotFoundError(f"Arquivo de dados não encontrado: {caminho}")

    if formato in (None, 'auto'):
        formato = detectar_formato(caminho)
        if formato is None:
            raise ValueError(
                f"Não foi possível detectar o formato de '{caminho}'. "
                f"Especifique 'formato' no YAML."
            )

    if formato not in FORMATOS:
        raise ValueError(
            f"Formato '{formato}' desconhecido. "
            f"Válidos: {sorted(FORMATOS)}"
        )

    return FORMATOS[formato](caminho, **kwargs)


def interpolar_em_grid(serie, t_grid):
    """Interpola uma série em um grid temporal comum. Ignora NaNs."""
    from scipy.interpolate import interp1d

    out = {}
    for chave in ('T_up_C', 'Z_int_m', 'O2_vol', 'CO2_vol', 'HRR_kW'):
        y = serie.get(chave)
        if y is None:
            out[chave] = None
            continue
        x = np.asarray(serie['t'], dtype=float)
        y = np.asarray(y, dtype=float)
        mask = np.isfinite(x) & np.isfinite(y)
        if mask.sum() < 2:
            out[chave] = np.full_like(t_grid, np.nan, dtype=float)
            continue
        f = interp1d(x[mask], y[mask], bounds_error=False, fill_value=np.nan)
        out[chave] = f(t_grid)
    return out
