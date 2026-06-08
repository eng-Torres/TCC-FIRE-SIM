"""
Banco de Dados de Materiais Combustíveis
Baseado em valores típicos de ensaios de cone calorimétrico (ISO 5660) e dados da literatura.
Referências:
- Schartel, B., Hull, T. R. (2007). Fire and Materials, 31, 327-354.
- SFPE Handbook, 5th Ed., Cap. 61 (Calorimetry).
- CFAST Technical Reference Guide (NIST SP 1026), Cap. 3 (Combustion).
"""

import numpy as np

# Banco de dados contendo:
# - 'tempo': array de tempo [s]
# - 'hrr': Heat Release Rate [kW] (para um objeto típico)
# - 'delta_Hc': calor de combustão efetivo [kJ/kg]
# - 'nu_o2': consumo estequiométrico de O2 [kg O2 / kg combustível]
# - 'y_co2': rendimento mássico de CO2 [kg CO2 / kg combustível]
# - 'y_co': rendimento mássico de CO [kg CO / kg combustível]
# Nota: os arrays de tempo/HRR foram gerados para representar um objeto inteiro,
# e não por unidade de área, para simplificar a integração inicial.
# Se desejar escalonar pela área, modifique a função de cálculo.

BANCO_MATERIAIS = {
    "espuma_poliuretano": {
        "tempo": np.array([0, 10, 20, 30, 40, 50, 60, 90, 120, 150, 180, 240, 300, 360, 420, 480, 540, 600]),
        "hrr":   np.array([0, 20, 80, 160, 250, 320, 380, 360, 310, 250, 190, 100, 50, 25, 10, 5, 2, 0]),
        "hrrpua_max": 800.0,  # kW/m²
        "delta_Hc": 20.0e3,   # kJ/kg (típico de PU flexível)
        "nu_o2": 2.2,         # kg O2 / kg fuel (base estequiométrica aproximada)
        "y_co2": 2.1,         # kg CO2 / kg fuel
        "y_co": 0.02          # kg CO / kg fuel (combustão bem ventilada)
    },
    "madeira": {
        "tempo": np.array([0, 30, 60, 90, 120, 180, 240, 300, 360, 420, 480, 540, 600]),
        "hrr":   np.array([0, 40, 100, 160, 200, 240, 230, 190, 140, 100, 60, 30, 10]),
        "hrrpua_max": 250.0,  # kW/m²
        "delta_Hc": 16.0e3,   # kJ/kg (madeira seca)
        "nu_o2": 1.4,         # kg O2 / kg fuel
        "y_co2": 1.5,         # kg CO2 / kg fuel
        "y_co": 0.005
    },
    "papelao": {
        "tempo": np.array([0, 10, 20, 30, 45, 60, 90, 120, 150, 180, 210, 240, 270, 300]),
        "hrr":   np.array([0, 10, 60, 140, 220, 280, 290, 240, 170, 100, 50, 20, 5, 0]),
        "hrrpua_max": 200.0,  # kW/m²
        "delta_Hc": 18.0e3,
        "nu_o2": 1.8,
        "y_co2": 1.9,
        "y_co": 0.01
    }
}

# Parâmetros padrão para uso com a curva t² quando o material não está no banco
ALFA_NFPA = {
    "lento": 0.00293,
    "medio": 0.01172,
    "rapido": 0.0469,
    "ultra_rapido": 0.1876
}

# Propriedades padrão de combustão para curva t² (material genérico)
COMBUSTAO_PADRAO = {
    "delta_Hc": 18.0e3,      # kJ/kg
    "nu_o2": 1.8,            # kg O2 / kg fuel
    "y_co2": 2.0,
    "y_co": 0.01
}

def obter_dados_material(nome_material: str) -> dict:
    """
    Retorna os dados de um material do banco de dados.
    Se não existir, retorna None (para usar curva t²).
    """
    return BANCO_MATERIAIS.get(nome_material, None)