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
        "delta_Hc": 20.0e3,   # kJ/kg (típico de PU flexível)
        "nu_o2": 2.2,         # kg O2 / kg fuel (base estequiométrica aproximada)
        "y_co2": 2.1,         # kg CO2 / kg fuel
        "y_co": 0.02          # kg CO / kg fuel (combustão bem ventilada)
    },
    "madeira": {
        "tempo": np.array([0, 30, 60, 90, 120, 180, 240, 300, 360, 420, 480, 540, 600]),
        "hrr":   np.array([0, 40, 100, 160, 200, 240, 230, 190, 140, 100, 60, 30, 10]),
        "delta_Hc": 16.0e3,   # kJ/kg (madeira seca)
        "nu_o2": 1.4,         # kg O2 / kg fuel
        "y_co2": 1.5,         # kg CO2 / kg fuel
        "y_co": 0.005
    },
    "papelao": {
        "tempo": np.array([0, 10, 20, 30, 45, 60, 90, 120, 150, 180, 210, 240, 270, 300]),
        "hrr":   np.array([0, 10, 60, 140, 220, 280, 290, 240, 170, 100, 50, 20, 5, 0]),
        "delta_Hc": 18.0e3,
        "nu_o2": 1.8,
        "y_co2": 1.9,
        "y_co": 0.01
    },
    "cesto_lixo": {
        "tempo": np.array([0, 200, 200.1, 600]),
        "hrr":   np.array([52.5, 52.5, 0, 0]),
        "delta_Hc": 20.0e3,       # kJ/kg (valor típico para papel/polietileno)
        "nu_o2": 2.0,             # kg O2 / kg fuel
        "y_co2": 2.0,             # kg CO2 / kg fuel
        "y_co": 0.02,             # kg CO / kg fuel
    },

    "armchair_h": {
        # Curva real do NIST "Arm Chair H" – substitui pelos dados do CSV se quiseres
        "tempo": np.array([0, 30, 60, 90, 120, 180, 240, 300, 360, 420, 480, 540, 600]),
        "hrr":   np.array([0, 100, 400, 800, 1200, 1600, 1733, 1500, 1200, 800, 400, 100, 0]),
        "delta_Hc": 20.0e3,
        "nu_o2": 2.2,
        "y_co2": 2.1,
        "y_co": 0.02,
    },
    "teste1_hrr_medida": {
    "tempo": np.array([
        0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 130, 140, 150, 160,
        170, 180, 190, 200, 210, 220, 230, 240, 250, 260, 270, 280, 290, 300, 310,
        320, 330, 340, 370, 380, 390, 400, 410, 420, 430, 440, 450,
        470, 480, 490, 500, 510, 520, 530, 540, 550, 560, 570, 580, 590, 600, 610,
        620, 630, 640, 650, 660, 670, 680, 690, 700, 710, 720, 730, 740, 750, 760,
        770, 780, 790, 800, 810, 820, 830, 840, 850, 860, 870, 880, 890, 900, 910,
        920, 930, 940, 950, 960, 970, 980, 990, 1000, 1010, 1020, 1030, 1040, 1050,
        1060, 1070, 1080, 1090, 1100, 1110, 1120, 1130, 1150, 1160, 1170, 1180,
        1200, 1210, 1220, 1230, 1240, 1250, 1260, 1270, 1280, 1290, 1300, 1310,
        1320, 1330, 1340, 1350, 1360, 1370, 1380, 1390, 1400, 1410, 1420, 1430, 1440,
        1450, 1460, 1470, 1480, 1490, 1500, 1510, 1520, 1530, 1540, 1550, 1560, 1570,
        1580, 1590, 1600, 1610, 1620, 1630, 1640, 1650, 1660, 1670, 1680, 1690, 1700,
        1710, 1720, 1730, 1740, 1750, 1760, 1770, 1780, 1800, 1810, 1820, 1830,
        1840, 1850, 1860, 1870, 1880, 1890, 1900, 1910, 1920, 1930, 1940, 1950, 1960,
        1970, 1980, 1990, 2000, 2010, 2030, 2040, 2050, 2060, 2070
    ]),
    "hrr": np.array([
        0, 0, 3.5416, 0, 0, 0, 0, 0, 1.861, 0, 0, 0, 8.5593, 0, 0, 0, 0, 5.2499,
        0.061592, 0, 0, 0, 21.965, 18.931, 0, 0, 0, 15.126, 40.292, 18.965, 45.799,
        49.978, 85.684, 151.78, 264.25, 402.93, 659.91, 1025.9, 1311.5, 1704.0,
        1862.4, 2050.3, 1996.3, 2342.0, 2433.7, 2372.8, 2174.9, 2031.5, 1980.1,
        1698.8, 1360.1, 1251.7, 1235.3, 1093.1, 943.45, 903.82, 767.71, 697.76,
        595.22, 555.45, 537.36, 474.03, 395.72, 392.37, 461.75, 523.8, 508.05,
        417.79, 391.04, 341.53, 320.23, 309.68, 307.03, 299.55, 288.01, 295.59,
        278.26, 272.72, 243.56, 222.45, 218.49, 219.8, 232.93, 223.51, 226.76,
        194.22, 176.28, 147.75, 154.77, 148.23, 172.85, 157.79, 156.74, 123.36,
        116.46, 110.94, 100.61, 133.96, 114.34, 123.17, 100.92, 90.888, 83.484,
        86.0, 106.75, 96.97, 108.91, 80.407, 67.142, 75.477, 58.296, 60.915,
        65.064, 87.976, 72.302, 68.97, 52.834, 62.319, 35.561, 47.763, 50.993,
        34.983, 61.423, 55.79, 61.738, 72.104, 50.996, 56.624, 39.403, 31.449,
        36.603, 24.812, 38.019, 59.145, 42.711, 57.239, 38.597, 37.081, 42.282,
        17.037, 30.815, 20.883, 8.501, 24.093, 9.1714, 31.765, 37.356, 25.591,
        40.321, 37.48, 20.951, 19.979, 13.765, 5.819, 19.945, 18.044, 15.646,
        36.45, 22.218, 27.923, 40.9, 24.751, 20.466, 21.752, 3.5021, 8.3958,
        19.446, 5.3156, 22.295, 28.699, 14.52, 28.545, 24.026, 4.4863, 17.953,
        0, 3.9531, 16.042, 9.4926, 22.071, 27.056, 11.27, 17.678, 27.048, 12.526,
        0, 13.666, 0, 0, 12.594, 14.107, 7.259, 14.466, 30.86, 18.922, 14.536,
        17.129, 18.18, 18.055, 0
    ]),

        
        #"tempo": np.array([0, 240, 300, 360, 420, 450, 480, 510, 540, 600, 660, 960, 1380, 2000]),
        #"hrr":   np.array([0,   0,  35, 360,1250,1400,1450,1400,1200, 800, 460, 120,   40,    0]),
        "delta_Hc": 20.0e3,       # valor médio (não crítico para a dinâmica térmica)
        "nu_o2": 2.0,             # valor típico
        "y_co2": 2.0,
        "y_co": 0.02,
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