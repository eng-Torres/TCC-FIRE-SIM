# src/edos_solver.py 
"""
Módulo de Dinâmica de Camadas
Resolve as Equações de Conservação de Massa, Energia e Espécies.
Baseado no SFPE Handbook (Cap. 36) e CFAST (NIST SP 1026).
Inclui correlação da pluma (Zukoski, McCaffrey ou Heskestad),
modelo de parede com condução transiente (N nós),
suporte a múltiplos objetos com ignição sequencial.
"""

import numpy as np
from scipy.integrate import solve_ivp
from src.fluidodinamica import calcular_fluxo_porta_duas_zonas
from src.combustao import calcular_hrr_total_e_propriedades, propriedades_combustao
from src.constantes import *                      # constantes físicas universais
from src.banco_materiais import MATERIAIS_CONSTRUCAO
from src.termica_parede import criar_parede, derivadas_parede

# ============================================================
# Parâmetros de simulação (ajustáveis)
# ============================================================
MATERIAL_PAREDE = 'gesso_16mm'   # chave em MATERIAIS_CONSTRUCAO
MAX_STEP = 0.1                        # passo máximo do integrador [s]
PONTOS_POR_SEGUNDO = 10               # densidade do t_eval
N_NOS = 12                             # número de nós na discretização da parede

# Propriedades da parede extraídas do banco de materiais
props_parede = MATERIAIS_CONSTRUCAO[MATERIAL_PAREDE]
K_W = props_parede['k']
RHO_W = props_parede['rho']
CP_W = props_parede['cp']
E_W = props_parede['espessura']
H_CONV = props_parede['h_conv']
epsilon = props_parede['emissividade']

# Inicializa o modelo de condução (parâmetros constantes da malha)
params_parede = criar_parede(
    N_nos=N_NOS,
    espessura_total=E_W,
    k=K_W,
    rho=RHO_W,
    cp=CP_W,
    h_conv_int=H_CONV,
    emissividade=epsilon,
    T_amb=T_AMB,
    h_conv_ext=0.005          # convecção exterior [kW/(m²·K)]
)

# -----------------------------------------------------------
# Correlações de pluma disponíveis
# -----------------------------------------------------------
def _pluma_zukoski(Q_c, z):
    """Zukoski (1981)"""
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0
    return 0.071 * (Q_c ** (1.0 / 3.0)) * (z ** (5.0 / 3.0))

def _pluma_heskestad(Q_dot_total, Q_c, z, D, T_inf=293.15, cp=1.005, rho_inf=1.2, g=9.81):
    """
    Heskestad plume entrainment – CFAST Eqs. (3.16)-(3.18).
    Usa Q_dot_total para a origem virtual e altura de chama,
    e Q_c para o termo convectivo da pluma.
    """
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0

    # Parâmetros adimensionais
    Q_star = Q_dot_total / (rho_inf * cp * T_inf * np.sqrt(g * D) * D**2)
    
    # Altura média da chama (Eq. 3.18)
    L = D * (-1.02 + 3.7 * Q_star**(2/5))
    L = max(0.01, L)   # não pode ser negativo ou zero

    # Origem virtual (Eq. 3.17) – note que usa Q_dot_total
    z0 = D * (-1.02 + 1.4 * Q_star**(2/5))

    # Constante comum
    const = 0.196 * (g * rho_inf**2 / (cp * T_inf))**(1/3)

    if z < L:
        # Abaixo da altura da chama: relação linear (Eq. 3.18)
        # Calcula o entrainment na altura L
        z_eff_L = max(L - z0, 0.01)
        m_e_L = const * Q_c**(1/3) * z_eff_L**(5/3)
        m_e_L *= 1.0 + (2.9 * Q_c**(2/3)) / ((np.sqrt(g) * cp * rho_inf * T_inf)**(2/3) * z_eff_L**(5/3))
        # Entrainment linear abaixo de L
        m_e = m_e_L * (z / L)
    else:
        # Acima da altura da chama: Eq. (3.16) completa
        z_eff = max(z - z0, 0.01)
        m_e = const * Q_c**(1/3) * z_eff**(5/3)
        m_e *= 1.0 + (2.9 * Q_c**(2/3)) / ((np.sqrt(g) * cp * rho_inf * T_inf)**(2/3) * z_eff**(5/3))

    return m_e

def _pluma_mccaffrey(Q_c, z):
    """McCaffrey (1979) – vazão mássica total [kg/s]."""
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0
    zeta = z / (Q_c ** 0.4)          # z / Q_c^{2/5}
    if zeta < 0.08:
        return 0.011 * Q_c * (zeta ** 0.566)
    elif zeta < 0.20:
        return 0.026 * Q_c * (zeta ** 0.909)
    else:
        return 0.124 * Q_c * (zeta ** 1.895)

def equacoes_sfpe(t, y, sala, aberturas, objetos_ativos, crescimento_t2, modelo_pluma='zukoski'):
    """
    Sistema de EDOs de conservação com múltiplos objetos.
    y = [m_u, H_u, m_o2, m_co2, m_co, T_w0, T_w1, ..., T_w(N-1)]
    modelo_pluma : 'zukoski', 'heskestad' ou 'mccaffrey'
    """
    # Extrair as variáveis 
    m_u   = y[0]
    H_u   = y[1]
    m_o2  = y[2]
    m_co2 = y[3]
    m_co  = y[4]
    T_parede = y[5:5+N_NOS]      # temperaturas dos nós da parede

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

    # Diâmetro fixo para a pluma de Heskestad (não afeta outras)
    D_fogo = 0.5

    # Escolha da função de pluma (usa o diâmetro fixo se for Heskestad)
    if modelo_pluma == 'zukoski':
        _calc_pluma = lambda Qc, z: _pluma_zukoski(Qc, z)
    elif modelo_pluma == 'mccaffrey':
        _calc_pluma = lambda Qc, z: _pluma_mccaffrey(Qc, z)
    else:
        _calc_pluma = lambda Qc, z: _pluma_heskestad(Q_dot, Q_c, z, D_fogo)

    Q_c = 0.7 * Q_dot
    m_dot_plume_total = _calc_pluma(Q_c, Z_int)
    m_dot_entrained = max(0.0, m_dot_plume_total - m_dot_fuel)

    # ------------------------------------------------------------
    # Ventilação por todas as aberturas (soma de fluxos)
    # ------------------------------------------------------------
    rho_u = P_ATM / (R_AR * T_u)
    m_dot_in = 0.0
    m_dot_out = 0.0
    for vent in aberturas:
        fluxo = calcular_fluxo_porta_duas_zonas(
            largura=vent["dimensoes"]["x"],
            altura_porta=vent["dimensoes"]["z"],
            altura_interface=Z_int,
            densidade_ar=RHO_AMB,
            densidade_fumaca=rho_u,
            C_d=0.68,
            m_dot_combustivel=m_dot_fuel
        )
        m_dot_in += fluxo["m_dot_in"]
        m_dot_out += fluxo["m_dot_out"]

    # Consumo e produção de espécies
    m_dot_o2_consumido = m_dot_fuel * nu_o2
    m_dot_co2_produzido = m_dot_fuel * y_co2
    m_dot_co_produzido  = m_dot_fuel * y_co

    # Limitação simples pelo O2 entranhado
    m_dot_o2_disp = m_dot_entrained * Y_O2_AMB
    if m_dot_o2_consumido > m_dot_o2_disp and m_dot_o2_disp > 0:
        fator = m_dot_o2_disp / m_dot_o2_consumido
        m_dot_fuel *= fator
        Q_dot *= fator
        Q_c = 0.7 * Q_dot
        # Recalcular a pluma e espécies
        m_dot_plume_total = _calc_pluma(Q_c, Z_int)
        m_dot_entrained = max(0.0, m_dot_plume_total - m_dot_fuel)
        m_dot_o2_consumido = m_dot_fuel * nu_o2
        m_dot_co2_produzido = m_dot_fuel * y_co2
        m_dot_co_produzido  = m_dot_fuel * y_co

    # Balanço das espécies na camada superior
    dm_o2_dt  = m_dot_entrained * Y_O2_AMB - m_dot_o2_consumido - m_dot_out * (m_o2 / m_u)
    dm_co2_dt = m_dot_entrained * Y_CO2_AMB + m_dot_co2_produzido - m_dot_out * (m_co2 / m_u)
    dm_co_dt  = m_dot_entrained * Y_CO_AMB  + m_dot_co_produzido  - m_dot_out * (m_co / m_u)

    # Evitar O2 negativo
    if m_o2 <= 0.0 and dm_o2_dt < 0.0:
        dm_o2_dt = 0.0

    # Balanço de massa total
    dm_dt = m_dot_plume_total - m_dot_out

    # Geometria para perda de calor
    L = sala["comprimento_x"]
    W = sala["profundidade_y"]
    H_sala = sala["altura_z"]
    A_teto = L * W
    A_paredes_sup = 2 * (L + W) * (H_sala - Z_int)
    A_contato = A_teto + A_paredes_sup

    # Perda de calor da camada superior para a superfície interna (nó 0)
    T_w0 = T_parede[0]
    Q_rad = epsilon * SIGMA * A_contato * (T_u**4 - T_w0**4)
    Q_conv = H_CONV * A_contato * (T_u - T_w0)
    Q_loss = Q_conv + Q_rad

    # Balanço de energia da camada superior
    dH_dt = (Q_c
             + m_dot_entrained * C_P * T_AMB
             - m_dot_out * C_P * T_u
             - Q_loss)

    # Derivadas da parede (chama o módulo termica_parede)
    dT_parede_dt = derivadas_parede(T_u, T_parede, params_parede)

    # As derivadas são: [massa, entalpia, O2, CO2, CO] + [N_NOS temperaturas]
    return [dm_dt, dH_dt, dm_o2_dt, dm_co2_dt, dm_co_dt] + list(dT_parede_dt)


def simular_incendio_sequencial(sala, aberturas,
                                objeto_inicial, outros_objetos,
                                objetos_geometria,
                                t_final=600, crescimento_t2="medio",
                                modelo_pluma='zukoski'):
    """
    Simulação completa com ignição sequencial.
    Retorna arrays de tempo e estados (t, y).
    modelo_pluma : 'zukoski', 'heskestad' ou 'mccaffrey'
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
    delta_T = 0.0
    T_u0 = T_AMB + delta_T
    rho_u0 = P_ATM / (R_AR * T_u0)
    m_u_0 = rho_u0 * V_u_0
    H_u_0 = m_u_0 * C_P * T_u0
    m_o2_0 = m_u_0 * Y_O2_AMB
    m_co2_0 = m_u_0 * Y_CO2_AMB
    m_co_0 = m_u_0 * Y_CO_AMB

    # Inicializar todas as temperaturas da parede a T_AMB
    T_parede0 = np.full(N_NOS, T_AMB)

    # Construir o vetor y0 
    y0 = np.concatenate([np.array([m_u_0, H_u_0, m_o2_0, m_co2_0, m_co_0]),
                         T_parede0])

    t_atual = 0.0
    y_atual = y0.copy()
    solucoes = []

    while t_atual < t_final and len(objetos_inativos) > 0:
        T_IGN_ATUAL = min(obj['T_ign'] for obj in objetos_inativos)

        def evento_ignicao(t, y, *args):
            m_u = y[0]
            H_u = y[1]
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
            
    # --- Continua a simulação até t_final se ainda houver tempo restante ---
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
    # ------------------------------------------------------------------------
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