"""
Módulo de Dinâmica de Camadas
Resolve as Equações de Conservação de Massa, Energia e Espécies.
Baseado no SFPE Handbook (Cap. 36) e CFAST (NIST SP 1026).
Inclui correlação da pluma (Zukoski ou Heskestad),
modelo de parede com inércia térmica e radiação,
suporte a múltiplos objetos com ignição sequencial,
e diâmetro da base da chama dinâmico e suave (HRRPUA + filtro).
"""

import numpy as np
from scipy.integrate import solve_ivp
from src.fluidodinamica import calcular_fluxo_porta_duas_zonas
from src.combustao import calcular_hrr_total_e_propriedades, propriedades_combustao
from src.constantes import *                      # constantes físicas universais
from src.banco_materiais import MATERIAIS_CONSTRUCAO

# ============================================================
# Parâmetros de simulação (ajustáveis)
# ============================================================
MATERIAL_PAREDE = 'concreto_normal'   # chave em MATERIAIS_CONSTRUCAO
TAU_D = 20.0                          # constante de tempo do filtro do diâmetro [s]
D_FOGO_INICIAL = 0.5                  # diâmetro inicial do fogo [m]
MAX_STEP = 0.1                        # passo máximo do integrador [s]
PONTOS_POR_SEGUNDO = 10               # densidade do t_eval

# Propriedades da parede extraídas do banco de materiais
props_parede = MATERIAIS_CONSTRUCAO[MATERIAL_PAREDE]
K_W = props_parede['k']
RHO_W = props_parede['rho']
CP_W = props_parede['cp']
E_W = props_parede['espessura']
H_CONV = props_parede['h_conv']
epsilon = props_parede['emissividade']

# Coeficiente auxiliar da parede (condução para o exterior)
C2 = K_W / (RHO_W * CP_W * E_W**2)     # [1/s]

# -----------------------------------------------------------
# Correlações de pluma disponíveis
# -----------------------------------------------------------
def _pluma_zukoski(Q_c, z):
    """Zukoski (1981)"""
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0
    return 0.071 * (Q_c ** (1.0 / 3.0)) * (z ** (5.0 / 3.0))

def _pluma_heskestad(Q_c, z, D):
    """Heskestad (1983) com transição suave na origem virtual."""
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0
    z0 = -1.02 * D + 0.083 * (Q_c ** 0.4)
    delta = max(0.1, abs(z0) * 0.1)
    if z <= z0 - delta:
        return 0.0
    elif z <= z0 + delta:
        t = (z - (z0 - delta)) / (2.0 * delta)
        f = t**2 * (3.0 - 2.0 * t)
    else:
        f = 1.0
    z_eff = max(z - z0, 0.01)
    m = 0.071 * (Q_c ** (1.0 / 3.0)) * (z_eff ** (5.0 / 3.0))
    m *= 1.0 + 0.026 * (Q_c ** (2.0 / 3.0)) / (z_eff ** (5.0 / 3.0))
    return f * m


def equacoes_sfpe(t, y, sala, aberturas, objetos_ativos, crescimento_t2, modelo_pluma='zukoski'):
    """
    Sistema de EDOs de conservação com múltiplos objetos.
    y = [m_u, H_u, m_o2, m_co2, m_co, T_w, D_fogo]
    modelo_pluma : 'zukoski' ou 'heskestad'
    """
    m_u, H_u, m_o2, m_co2, m_co, T_w, D_fogo = y

    m_u = max(m_u, 0.001)
    m_o2 = max(m_o2, 0.0)
    m_co2 = max(m_co2, 0.0)
    m_co = max(m_co, 0.0)

    # Propriedades da camada superior
    T_u = H_u / (m_u * C_P)
    T_u = max(T_u, 200.0)                         # proteção contra temperaturas irreais
    V_u = (m_u * R_AR * T_u) / P_ATM
    area_sala = sala["area_piso"]
    Z_int = sala["altura_z"] - (V_u / area_sala)
    Z_int = max(Z_int, 0.0)

    # Obter HRR total e propriedades ponderadas
    Q_dot, m_dot_fuel, delta_Hc, nu_o2, y_co2, y_co = \
        calcular_hrr_total_e_propriedades(t, objetos_ativos)

    # Cálculo do diâmetro instantâneo com base na HRR e HRRPUA
    total_area_objetos = sum(obj.get('area', 0.0) for obj in objetos_ativos)
    if total_area_objetos <= 0.0:
        D_calc = 0.5
    else:
        hrrpua_ponderado = 0.0
        for obj in objetos_ativos:
            props = propriedades_combustao(obj['nome_material'])
            hrrpua_ponderado += obj.get('area', 0.0) * props.get('hrrpua_max', 500.0)
        hrrpua_ponderado /= total_area_objetos

        area_ativa = Q_dot / max(hrrpua_ponderado, 1.0)
        area_ativa = min(area_ativa, total_area_objetos)
        area_ativa = max(area_ativa, 0.1)
        D_calc = np.sqrt(4.0 * area_ativa / np.pi)

    # Filtro passa‑baixo para o diâmetro
    dD_dt = (D_calc - D_fogo) / TAU_D

    # Escolha da função de pluma (usa o diâmetro filtrado)
    if modelo_pluma == 'zukoski':
        _calc_pluma = lambda Qc, z: _pluma_zukoski(Qc, z)
    else:
        _calc_pluma = lambda Qc, z: _pluma_heskestad(Qc, z, D_fogo)

    Q_c = 0.7 * Q_dot
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

    # Perda de calor da camada superior para a parede (convecção + radiação)
    Q_rad = epsilon * SIGMA * A_contato * (T_u**4 - T_w**4)
    Q_conv = H_CONV * A_contato * (T_u - T_w)
    Q_loss = Q_conv + Q_rad
    q_dot_w = H_CONV * (T_u - T_w) + epsilon * SIGMA * (T_u**4 - T_w**4)

    # Balanço de energia da camada superior
    dH_dt = (Q_c
             + m_dot_entrained * C_P * T_AMB
             - m_dot_out * C_P * T_u
             - Q_loss)

    # Balanço da parede
    C1 = q_dot_w / (RHO_W * CP_W * E_W)          # ganho de calor [1/s]
    dTw_dt = C1 - C2 * (T_w - T_AMB)

    return [dm_dt, dH_dt, dm_o2_dt, dm_co2_dt, dm_co_dt, dTw_dt, dD_dt]


def simular_incendio_sequencial(sala, aberturas,
                                objeto_inicial, outros_objetos,
                                objetos_geometria,
                                t_final=600, crescimento_t2="medio",
                                modelo_pluma='zukoski'):
    """
    Simulação completa com ignição sequencial.
    Retorna arrays de tempo e estados (t, y).
    modelo_pluma : 'zukoski' ou 'heskestad'
    """
    # Dicionário para lookup da área de cada material
    area_por_material = {obj['nome_material']: obj['area'] for obj in objetos_geometria}
    
    nome_inicial = objeto_inicial['nome_material']
    area_inicial = area_por_material.get(nome_inicial, 2.8)   # fallback (cama padrão)
    objetos_ativos = [{'nome_material': nome_inicial, 't_ign': 0.0, 'area': area_inicial}]

    objetos_inativos = []
    for mat in outros_objetos:
        nome = mat['nome_material']
        area = area_por_material.get(nome, 1.0)
        objetos_inativos.append({
            'nome_material': nome,
            't_ign': None,
            'T_ign': mat['T_ign'],
            'area': area
        })

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
    D_fogo0 = D_FOGO_INICIAL

    y0 = np.array([m_u_0, H_u_0, m_o2_0, m_co2_0, m_co_0, T_w0, D_fogo0])

    t_atual = 0.0
    y_atual = y0.copy()
    solucoes = []

    while t_atual < t_final and len(objetos_inativos) > 0:
        T_IGN_ATUAL = min(obj['T_ign'] for obj in objetos_inativos)

        def evento_ignicao(t, y, *args):
            m_u, H_u = y[0], y[1]
            T_u = H_u / (m_u * C_P)
            return T_u - T_IGN_ATUAL
        evento_ignicao.terminal = True
        evento_ignicao.direction = 1

        t_span = (t_atual, t_final)
        num_points = max(10, int(t_span[1] - t_span[0]) * PONTOS_POR_SEGUNDO)
        t_eval = np.linspace(t_span[0], t_span[1], num_points)

        sol = solve_ivp(
            equacoes_sfpe,
            t_span,
            y_atual,
            args=(sala, aberturas, objetos_ativos, crescimento_t2, modelo_pluma),
            events=evento_ignicao,
            method='Radau',
            max_step=MAX_STEP,
            t_eval=t_eval
        )
        solucoes.append(sol)
        t_atual = sol.t[-1]
        y_atual = sol.y[:, -1]

        if sol.t_events[0].size > 0 and len(objetos_inativos) > 0:
            obj_ignitado = min(objetos_inativos, key=lambda x: x['T_ign'])
            objetos_inativos.remove(obj_ignitado)
            obj_ignitado['t_ign'] = t_atual
            objetos_ativos.append(obj_ignitado)
            print(f"Ignição de {obj_ignitado['nome_material']} em t = {t_atual:.2f} s")
        else:
            if t_atual < t_final:
                t_span_rest = (t_atual, t_final)
                num_points = max(10, int(t_span_rest[1] - t_span_rest[0]) * PONTOS_POR_SEGUNDO)
                t_eval_rest = np.linspace(t_span_rest[0], t_span_rest[1], num_points)
                sol_rest = solve_ivp(
                    equacoes_sfpe,
                    t_span_rest,
                    y_atual,
                    args=(sala, aberturas, objetos_ativos, crescimento_t2, modelo_pluma),
                    method='Radau',
                    max_step=MAX_STEP,
                    t_eval=t_eval_rest
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