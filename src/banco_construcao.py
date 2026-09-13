"""
Banco de propriedades de materiais de construção (paredes, teto, piso).
"""

MATERIAIS_CONSTRUCAO = {
    "concreto_normal": {
        "k": 1.75e-3,
        "rho": 2200.0,
        "cp": 1.0,
        "espessura": 0.15,
        "h_conv": 0.015,
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
        "k": 45.0e-3,
        "rho": 7850.0,
        "cp": 0.5,
        "espessura": 0.005,
        "h_conv": 0.02,
        "emissividade": 0.3
    },
    "gesso_16mm": {          # Babrauskas (1984) — 16 mm Type X
        "k": 2.0e-4,         # kW/(m·K)
        "rho": 800.0,        # kg/m³
        "cp": 1.0,           # kJ/(kg·K)
        "espessura": 0.016,  # 16 mm
        "h_conv": 0.008,     # kW/(m²·K) — dentro da faixa CFAST
        "emissividade": 0.9  # gesso é altamente emissivo
    }
}