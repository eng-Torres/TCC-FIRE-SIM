"""
Módulo de Dinâmica de Camadas
Resolve as Equações de Conservação de Massa e Energia (SFPE Handbook).
Inclui modelo não-adiabático e escoamento bidirecional.
"""
import numpy as np
from scipy.integrate import solve_ivp
from tcc.src.fluidodinamica import calcular_fluxo_porta_duas_zonas
from tcc.src.combustão import calcular_hrr

# Constantes Termodinâmicas
C_P = 1.005      # Calor específico a pressão constante (kJ/kg.K)
R_AR = 0.287     # Constante do gás para o ar (kJ/kg.K)
P_ATM = 101.325  # Pressão atmosférica (kPa)
T_AMB = 293.15   # Temperatura ambiente (K) - equivale a 20°C
RHO_AMB = 1.2    # Densidade do ar ambiente (kg/m3)

def equacoes_sfpe(t, y, sala, aberturas, crescimento_fogo):
    """
    Sistema de EDOs de conservação.
    y[0] = m_u (Massa da camada superior em kg)
    y[1] = H_u (Entalpia total da camada superior em kJ)
    """
    m_u, H_u = y
    
    # Prevenção numérico-matemática para evitar singularidades
    m_u = max(m_u, 0.001)
    
    # 1. Recuperar propriedades a partir da massa e entalpia
    T_u = H_u / (m_u * C_P)
    
    # Lei dos Gases Ideais para achar o volume da fumaça
    V_u = (m_u * R_AR * T_u) / P_ATM
    
    # Altura da interface da fumaça (Z_int). Começa no teto e desce.
    area_sala = sala["area_piso"]
    Z_int = sala["altura_z"] - (V_u / area_sala)
    Z_int = max(Z_int, 0.0) 
    
    # 2. Calcular a Potência do Fogo (HRR) e a Fração Convectiva
    Q_dot = calcular_hrr(t, crescimento=crescimento_fogo, limite_max_kw=3000)
    Q_c = 0.7 * Q_dot 
    
    # 3. Calcular Fluxo de Massa da Pluma (Equação de Zukoski)
    if Z_int > 0.01:
        m_dot_p = 0.071 * (Q_c ** (1/3)) * (Z_int ** (5/3))
    else:
        m_dot_p = 0.0
        
    # 4. Calcular Ventilação Bidirecional da Porta
    # A função calcular_fluxo_porta_duas_zonas deve retornar entrada, saída e plano neutro
    porta = aberturas[0]
    m_dot_in, m_dot_out, Z_n = calcular_fluxo_porta_duas_zonas(
        largura=porta["dimensoes"]["x"],
        altura_porta=porta["dimensoes"]["z"],
        altura_interface=Z_int,
        densidade_ar=RHO_AMB,
        densidade_fumaca=rho_u,
        C_d=0.68,
        m_dot_combustivel=0.0                                                                  
    )

    # 5. Calcular Perda de Calor para as Superfícies (Modelo Não-Adiabático)
    L = sala["dimensoes"]["x"]
    W = sala["dimensoes"]["y"]
    H_sala = sala["altura_z"]
    
    A_teto = L * W
    A_paredes_sup = 2 * (L + W) * (H_sala - Z_int)
    A_contato = A_teto + A_paredes_sup
    
    h_eff = 0.03 # Coeficiente global de transferência térmica (kW/m2.K)
    Q_loss = h_eff * A_contato * (T_u - T_AMB)
        
    # 6. As Equações Diferenciais de Conservação
    dm_dt = m_dot_p - m_dot_out
    dH_dt = Q_c + (m_dot_p * C_P * T_AMB) - (m_dot_out * C_P * T_u) - Q_loss
    
    return [dm_dt, dH_dt]

def simular_incendio(sala, aberturas, crescimento="medio", t_final=600):
    """
    Executa o integrador numérico SciPy para resolver a evolução do incêndio.
    """
    # Condições Iniciais (Pequeno volume inicial a 20°C)
    V_u_0 = 0.01  # m³
    delta_T_inicial = 0.1  # pequeno superaquecimento para evitar densidade igual
    T_u0 = T_AMB + delta_T_inicial
    rho_u0 = P_ATM / (R_AR * T_u0)   # ou use rho_amb * T_amb / T_u0
    m_u_0 = rho_u0 * V_u_0
    H_u_0 = m_u_0 * C_P * T_u0
    y0 = [m_u_0, H_u_0]
    t_span = (0, t_final)
    t_eval = np.linspace(0, t_final, t_final)
    
    sol = solve_ivp(
        equacoes_sfpe, 
        t_span, 
        y0, 
        t_eval=t_eval, 
        args=(sala, aberturas, crescimento),
        method='Radau'
    )
    
    return sol