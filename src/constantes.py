"""
Constantes físicas universais e condições ambiente padrão.
Valores comuns a todos os módulos do simulador.
"""

# Termodinâmicas do ar
C_P = 1.005          # kJ/(kg·K) – calor específico a pressão constante
R_AR = 0.287         # kJ/(kg·K) – constante do gás para o ar
P_ATM = 101.325      # kPa – pressão atmosférica padrão

# Condições ambiente de referência
T_AMB = 293.15       # K (20 °C)
RHO_AMB = 1.2        # kg/m³ – densidade do ar ambiente

# Frações mássicas no ar ambiente
Y_O2_AMB = 0.23
Y_CO2_AMB = 0.0
Y_CO_AMB = 0.0

SIGMA = 5.67e-11   # kW/(m²·K⁴) – constante de Stefan-Boltzmann