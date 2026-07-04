"""
Banco de propriedades de materiais de construção (paredes, teto, piso).
Cada entrada é um dicionário com as propriedades térmicas necessárias
para o modelo de parede com inércia (lumped).
"""

MATERIAIS_CONSTRUCAO = {
    "concreto_normal": {
        "k": 1.75e-3,        # kW/(m·K)
        "rho": 2200.0,       # kg/m³
        "cp": 1.0,           # kJ/(kg·K)
        "espessura": 0.15,   # m
        "h_conv": 0.015,     # kW/(m²·K)
        "emissividade": 0.9
    },
    "gesso_cartonado": {
        "k": 2.0e-4,
        "rho": 800.0,
        "cp": 1.0,
        "espessura": 0.0125,
        "h_conv": 0.015,
        "emissividade": 0.9
    },
    "aco": {
        "k": 45.0e-3,        # kW/(m·K)  (45 W/mK)
        "rho": 7850.0,
        "cp": 0.5,
        "espessura": 0.005,
        "h_conv": 0.02,
        "emissividade": 0.3
    },
    "gesso_16mm": {          # ou "gesso_cartonado" se preferires
        "k": 2.0e-4,         # mantém
        "rho": 800.0,        # mantém
        "cp": 1.0,           # mantém
        "espessura": 0.0125, # 12.5 mm (padrão)
        "h_conv": 0.008,     # reduzido de 0.015
        "emissividade": 0.5  # reduzida de 0.9
    }
}