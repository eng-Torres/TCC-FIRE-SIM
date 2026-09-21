"""
Módulo de Física do Simulador de Incêndio
Cálculo da Taxa de Liberação de Calor (HRR) e propriedades de combustão.
Suporta curvas de materiais reais (cone calorimétrico) e o modelo t² da NFPA.
"""

import numpy as np
from src.banco_materiais import obter_dados_material, ALFA_NFPA, COMBUSTAO_PADRAO


def calcular_hrr(tempo_s: float, nome_material: str = "generico",
                 crescimento_t2: str = "medio") -> float:
    """Calcula a HRR [kW] em um instante t para um único material."""
    dados = obter_dados_material(nome_material)
    if dados is not None:
        tempos = dados["tempo"]
        hrrs = dados["hrr"]
        return float(np.interp(tempo_s, tempos, hrrs, left=0.0, right=0.0))
    else:
        alfa = ALFA_NFPA.get(crescimento_t2, ALFA_NFPA["medio"])
        return alfa * (tempo_s ** 2)


def propriedades_combustao(nome_material: str) -> dict:
    dados = obter_dados_material(nome_material)
    if dados is not None:
        return {
            "delta_Hc": dados["delta_Hc"],
            "nu_o2":    dados["nu_o2"],
            "y_co2":    dados["y_co2"],
            "y_co":     dados["y_co"],
        }
    else:
        return COMBUSTAO_PADRAO.copy()

HRR_MIN_PLUME = 2.50   # kW — abaixo disso, a pluma não é acionada
                       # justificativa: correlações de Heskestad/Zukoski/McCaffrey
                       # foram ajustadas para chamas com HRR >> 10 kW.
                       # Abaixo desse valor, a pluma não é turbulenta e a
                       # correlação perde validade física.

def calcular_hrr_total_e_propriedades(t_global, objetos_ativos):
    """
    Calcula HRR total, vazão mássica total, propriedades ponderadas
    e o diâmetro equivalente da fonte de fogo.

    Cada objeto ativo é um dict com:
        'nome_material': str
        't_ign': float
        'area':  float (m², opcional; default 0.25 m²)

    Retorna:
        hrr_total, m_fuel_total, delta_Hc_pond, nu_o2_pond,
        y_co2_pond, y_co_pond, D_fogo
    """
    hrr_total    = 0.0
    m_fuel_total = 0.0
    delta_Hc_num = 0.0
    nu_o2_num    = 0.0
    y_co2_num    = 0.0
    y_co_num     = 0.0
    area_hrr_num = 0.0   # Σ HRR_i · A_i  (para ponderar a área pela HRR)

    for obj in objetos_ativos:
        t_local = t_global - obj['t_ign']
        if t_local < 0:
            continue

        props   = propriedades_combustao(obj['nome_material'])
        hrr_obj = calcular_hrr(t_local, nome_material=obj['nome_material'])
        area_obj = float(obj.get('area', 0.25))

        m_fuel_obj = hrr_obj / props['delta_Hc'] if props['delta_Hc'] > 0 else 0.0

        hrr_total    += hrr_obj
        m_fuel_total += m_fuel_obj
        delta_Hc_num += m_fuel_obj * props['delta_Hc']
        nu_o2_num    += m_fuel_obj * props['nu_o2']
        y_co2_num    += m_fuel_obj * props['y_co2']
        y_co_num     += m_fuel_obj * props['y_co']
        area_hrr_num += hrr_obj * area_obj

    # ---- Limiar de pluma ----
    if hrr_total < HRR_MIN_PLUME:
        hrr_total    = 0.0
        m_fuel_total = 0.0
        delta_Hc_num = 0.0
        nu_o2_num    = 0.0
        y_co2_num    = 0.0
        y_co_num     = 0.0
        area_hrr_num = 0.0

    if m_fuel_total > 0:
        delta_Hc_pond = delta_Hc_num / m_fuel_total
        nu_o2_pond    = nu_o2_num    / m_fuel_total
        y_co2_pond    = y_co2_num    / m_fuel_total
        y_co_pond     = y_co_num     / m_fuel_total
    else:
        p = COMBUSTAO_PADRAO
        delta_Hc_pond = p['delta_Hc']
        nu_o2_pond    = p['nu_o2']
        y_co2_pond    = p['y_co2']
        y_co_pond     = p['y_co']

    # ---- Diâmetro equivalente da fonte (círculo de mesma área efetiva) ----
    if hrr_total > 0:
        A_ef = area_hrr_num / hrr_total
    else:
        A_ef = np.pi * (0.5 ** 2) / 4.0   
    D_fogo = max(0.10, float(np.sqrt(4.0 * A_ef / np.pi)))

    return (hrr_total, m_fuel_total, delta_Hc_pond,
            nu_o2_pond, y_co2_pond, y_co_pond, D_fogo)