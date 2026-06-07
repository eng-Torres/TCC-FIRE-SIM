"""
Módulo de Dinâmica de Camadas
Resolve as Equações de Conservação de Massa, Energia e Espécies.
Baseado no SFPE Handbook (Cap. 36) e CFAST (NIST SP 1026).
Inclui correlação da pluma (Zukoski ou Heskestad),
modelo de parede com inércia térmica,
suporte a múltiplos objetos com ignição sequencial.
"""

import numpy as np
from scipy.integrate import solve_ivp
from src.fluidodinamica import calcular_fluxo_porta_duas_zonas
from src.combustao import calcular_hrr_total_e_propriedades

# Constantes Termodinâmicas
C_P = 1.005          # kJ/(kg·K)
R_AR = 0.287         # kJ/(kg·K)
P_ATM = 101.325      # kPa
T_AMB = 293.15       # K (20 °C)
RHO_AMB = 1.2        # kg/m³
Y_O2_AMB = 0.23      # fração mássica de O2 no ar ambiente
Y_CO2_AMB = 0.0
Y_CO_AMB = 0.0

# Propriedades da parede (concreto normal)
K_W = 1.75e-3         # kW/(m·K)   (1.75 W/(m·K))
RHO_W = 2200.0        # kg/m³
CP_W = 1.0            # kJ/(kg·K)
E_W = 0.15            # m (espessura)
H_CONV = 0.015         # kW/(m²·K)  (coeficiente de convecção natural)



# -----------------------------------------------------------
# Correlações de pluma disponíveis
# -----------------------------------------------------------
def _pluma_zukoski(Q_c, z):
    """Zukoski (1981)"""
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0
    return 0.071 * (Q_c ** (1.0 / 3.0)) * (z ** (5.0 / 3.0))

def _pluma_heskestad(Q_c, z, D=1.89):
    """Heskestad (1983)"""
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0
    z0 = -1.02 * D + 0.083 * (Q_c ** 0.4)
    if z <= z0:
        return 0.0
    m = 0.071 * (Q_c ** (1.0 / 3.0)) * ((z - z0) ** (5.0 / 3.0))
    m *= 1.0 + 0.026 * (Q_c ** (2.0 / 3.0)) / ((z - z0) ** (5.0 / 3.0))
    return m


def equacoes_sfpe(t, y, sala, aberturas, objetos_ativos, crescimento_t2, modelo_pluma='zukoski'):
    """
    Sistema de EDOs de conservação com múltiplos objetos.
    y = [m_u, H_u, m_o2, m_co2, m_co, T_w]
    modelo_pluma : 'zukoski' ou 'heskestad'
    """
    m_u, H_u, m_o2, m_co2, m_co, T_w = y

    m_u = max(m_u, 0.001)
    m_o2 = max(m_o2, 0.0)
    m_co2 = max(m_co2, 0.0)
    m_co = max(m_co, 0.0)

    # Propriedades da camada superior
    T_u = H_u / (m_u * C_P)
    V_u = (m_u * R_AR * T_u) / P_ATM
    area_sala = sala["area_piso"]
    Z_int = sala["altura_z"] - (V_u / area_sala)
    Z_int = max(Z_int, 0.0)

    # Obter HRR total e propriedades ponderadas
    Q_dot, m_dot_fuel, delta_Hc, nu_o2, y_co2, y_co = \
        calcular_hrr_total_e_propriedades(t, objetos_ativos)

    # Escolha da função de pluma
    if modelo_pluma == 'zukoski':
        _calc_pluma = lambda Qc, z: _pluma_zukoski(Qc, z)
    else:
        _calc_pluma = lambda Qc, z: _pluma_heskestad(Qc, z)

    Q_c = 0.7 * Q_dot
    # Pluma – vazão total
    m_dot_plume_total = _calc_pluma(Q_c, Z_int)
    m_dot_entrained = max(0.0, m_dot_plume_total - m_dot_fuel)

    # Ventilação pela porta
    porta = aberturas[0]
    rho_u = P_ATM / (R_AR * T_u)
    fluxos_porta = calcular_fluxo_porta_duas_zonas(
        largura=porta["dimensoes"]["x"],
        altura_porta=porta["dimensoes"]["z"],
        altura_interface=Z_int,
        densidade_ar=RHO_AMB,
        densidade_fumaca=rho_u,
        C_d=0.68,
        m_dot_combustivel=m_dot_fuel
    )
    m_dot_in = fluxos_porta["m_dot_in"]
    m_dot_out = fluxos_porta["m_dot_out"]
    Z_n = fluxos_porta["z_n"]

    # Consumo e produção de espécies
    m_dot_o2_consumido = m_dot_fuel * nu_o2
    m_dot_co2_produzido = m_dot_fuel * y_co2
    m_dot_co_produzido  = m_dot_fuel * y_co

    # Limitação por oxigénio na pluma
    m_dot_o2_disp = m_dot_entrained * Y_O2_AMB
    if m_dot_o2_consumido > m_dot_o2_disp and m_dot_o2_disp > 0:
        fator = m_dot_o2_disp / m_dot_o2_consumido
        m_dot_fuel *= fator
        Q_dot *= fator
        Q_c = 0.7 * Q_dot
        m_dot_plume_total = _calc_pluma(Q_c, Z_int)
        m_dot_entrained = max(0.0, m_dot_plume_total - m_dot_fuel)
        m_dot_o2_consumido = m_dot_fuel * nu_o2
        m_dot_co2_produzido = m_dot_fuel * y_co2
        m_dot_co_produzido  = m_dot_fuel * y_co

    # -----------------------------------------------------------
    # Extinção suave
    # -----------------------------------------------------------
    Qc_lim_inf = 0.5      # kW
    Qc_lim_sup = 5.0      # kW
    deltaT_lim = 3.0      # K

    if t > 1.0:
        fator_Qc = np.clip((Q_c - Qc_lim_inf) / (Qc_lim_sup - Qc_lim_inf), 0.0, 1.0)
        fator_dT = np.clip((T_u - T_AMB) / deltaT_lim, 0.0, 1.0)
        fator_extincao = min(fator_Qc, fator_dT)

        Q_c *= fator_extincao
        Q_dot *= fator_extincao

        m_dot_plume_total = _calc_pluma(Q_c, Z_int)

        if delta_Hc > 0:
            m_dot_fuel = Q_dot / delta_Hc
        else:
            m_dot_fuel = 0.0

        m_dot_entrained = max(0.0, m_dot_plume_total - m_dot_fuel)
        m_dot_o2_consumido = m_dot_fuel * nu_o2
        m_dot_co2_produzido = m_dot_fuel * y_co2
        m_dot_co_produzido  = m_dot_fuel * y_co
    # -----------------------------------------------------------

    # Balanço das espécies na camada superior
    dm_o2_dt  = m_dot_entrained * Y_O2_AMB - m_dot_o2_consumido - m_dot_out * (m_o2 / m_u)
    dm_co2_dt = m_dot_entrained * Y_CO2_AMB + m_dot_co2_produzido - m_dot_out * (m_co2 / m_u)
    dm_co_dt  = m_dot_entrained * Y_CO_AMB  + m_dot_co_produzido  - m_dot_out * (m_co / m_u)

    # Balanço de massa total
    dm_dt = m_dot_plume_total - m_dot_out

    # Geometria para perda de calor
    L = sala["comprimento_x"]
    W = sala["profundidade_y"]
    H_sala = sala["altura_z"]
    A_teto = L * W
    A_paredes_sup = 2 * (L + W) * (H_sala - Z_int)
    A_contato = A_teto + A_paredes_sup

    # Perda de calor da camada superior para a parede
    epsilon = 0.9
    sigma = 5.67e-11   # kW/(m²·K⁴)
    Q_rad = epsilon * sigma * A_contato * (T_u**4 - T_w**4)
    Q_conv = H_CONV * A_contato * (T_u - T_w)
    Q_loss = Q_conv + Q_rad  # kW
    q_dot_w = H_CONV * (T_u - T_w) + epsilon * sigma * (T_u**4 - T_w**4)   # kW/m²


    # Balanço de energia da camada superior
    dH_dt = (Q_c
             + m_dot_entrained * C_P * T_AMB
             - m_dot_out * C_P * T_u
             - Q_loss)

    # Balanço da parede (modelo lumped)
    # Coeficientes da equação da parede (já divididos pela capacitância)
    C1 = q_dot_w / (RHO_W * CP_W * E_W)       # 1/s
    C2 = K_W / (RHO_W * CP_W * E_W**2)        # 1/s
    # Fluxo total para a parede (convecção + radiação)
    dTw_dt = C1 - C2 * (T_w - T_AMB)

    return [dm_dt, dH_dt, dm_o2_dt, dm_co2_dt, dm_co_dt, dTw_dt]


def simular_incendio_sequencial(sala, aberturas,
                                objeto_inicial, outros_objetos,
                                t_final=600, crescimento_t2="medio",
                                modelo_pluma='zukoski'):
    """
    Simulação completa com ignição sequencial.
    Retorna arrays de tempo e estados (t, y).
    modelo_pluma : 'zukoski' ou 'heskestad'
    """
    objetos_ativos = [{'nome_material': objeto_inicial['nome_material'], 't_ign': 0.0}]
    objetos_inativos = [{'nome_material': mat, 't_ign': None} for mat in outros_objetos]

    # Condições iniciais
    V_u_0 = 0.01
    delta_T = 0.1
    T_u0 = T_AMB + delta_T
    rho_u0 = P_ATM / (R_AR * T_u0)
    m_u_0 = rho_u0 * V_u_0
    H_u_0 = m_u_0 * C_P * T_u0
    m_o2_0 = m_u_0 * Y_O2_AMB
    m_co2_0 = m_u_0 * Y_CO2_AMB
    m_co_0 = m_u_0 * Y_CO_AMB
    T_w0 = T_AMB

    y0 = np.array([m_u_0, H_u_0, m_o2_0, m_co2_0, m_co_0, T_w0])

    t_atual = 0.0
    y_atual = y0.copy()
    solucoes = []

    T_IGN = 573.15  # 300 °C

    while t_atual < t_final and len(objetos_inativos) > 0:
        def evento_ignicao(t, y, *args):
            m_u, H_u = y[0], y[1]
            T_u = H_u / (m_u * C_P)
            return T_u - T_IGN
        evento_ignicao.terminal = True
        evento_ignicao.direction = 1

        t_span = (t_atual, t_final)
        sol = solve_ivp(
            equacoes_sfpe,
            t_span,
            y_atual,
            args=(sala, aberturas, objetos_ativos, crescimento_t2, modelo_pluma),
            events=evento_ignicao,
            method='Radau',
            max_step=1.0
        )
        solucoes.append(sol)
        t_atual = sol.t[-1]
        y_atual = sol.y[:, -1]

        if sol.t_events[0].size > 0 and len(objetos_inativos) > 0:
            novo_obj = objetos_inativos.pop(0)
            novo_obj['t_ign'] = t_atual
            objetos_ativos.append(novo_obj)
            print(f"Ignição de {novo_obj['nome_material']} em t = {t_atual:.2f} s")
        else:
            if t_atual < t_final:
                t_span_rest = (t_atual, t_final)
                sol_rest = solve_ivp(
                    equacoes_sfpe,
                    t_span_rest,
                    y_atual,
                    args=(sala, aberturas, objetos_ativos, crescimento_t2, modelo_pluma),
                    method='Radau',
                    max_step=1.0
                )
                solucoes.append(sol_rest)
                t_atual = t_final
            break

    # Concatenação segura do histórico
    t_hist = [np.array([0.0])]
    y_hist = [y0.reshape(-1, 1)]
    for sol_i in solucoes:
        if sol_i.t.size > 1:
            t_hist.append(sol_i.t[1:])
            y_hist.append(sol_i.y[:, 1:])
    t_total = np.concatenate(t_hist)
    y_total = np.concatenate(y_hist, axis=1)
    return t_total, y_total