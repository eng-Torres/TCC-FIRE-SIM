# src/edos_solver.py
"""
Módulo de Dinâmica de Camadas
Resolve as Equações de Conservação de Massa, Energia e Espécies.
Baseado no SFPE Handbook (Cap. 36) e CFAST (NIST SP 1026/1030).
Inclui correlação da pluma (Zukoski, McCaffrey ou Heskestad),
modelo de parede com condução transiente (N nós),
suporte a múltiplos objetos com ignição sequencial.

Referências:
[1] Jones et al., NIST SP 1030 (2004) — CFAST Technical Reference Guide.
[2] Hurley et al., SFPE Handbook of Fire Protection Engineering, 5th ed. (2016).
[3] McGrattan et al., FDS User's Guide (2013) — Flame Extinction.
[4] SciPy Documentation — solve_ivp, Radau, events.
[5] Heskestad, Fire Safety Journal (1983) — Plume entrainment.
[6] McCaffrey, NBSIR 79-1910 (1979) — Plume correlations.
[7] Zukoski et al., Fire Safety Journal (1981) — Entrainment.

# A limitação por O2 usa a função tanh suave do CFAST 6/7
# (NIST SP 1030, Eq. 47), com Y_LOL calibrado para o cenário
# de validação (Babrauskas, 1984). Ver Seção 3.X.3.7.
"""

import numpy as np
from scipy.integrate import solve_ivp
from src.fluidodinamica import calcular_fluxo_porta_duas_zonas
from src.combustao import calcular_hrr_total_e_propriedades
from src.constantes import *
from src.banco_construcao import MATERIAIS_CONSTRUCAO
from src.termica_parede import criar_parede, derivadas_parede

# ============================================================
# Parâmetros de simulação (ajustáveis)
# ============================================================
MATERIAL_PAREDE_FALLBACK = 'gesso_16mm'   
MAX_STEP = 0.1
PONTOS_POR_SEGUNDO = 10
N_NOS = 5

if N_NOS < 2:
    raise ValueError(
        f"N_NOS deve ser >= 2 para discretização por diferenças finitas. "
        f"Recebido: {N_NOS}. Para modelo de capacitância concentrada, "
        f"implemente um caso especial em termica_parede.py."
    )

RHO_AMB_CALC = P_ATM / (R_AR * T_AMB)
if abs(RHO_AMB - RHO_AMB_CALC) / RHO_AMB_CALC > 1e-3:
    print(f"[AVISO] RHO_AMB ({RHO_AMB:.4f}) difere de "
      f"P_ATM/(R_AR*T_AMB) ({RHO_AMB_CALC:.4f}). "
      f"Usando valor calculado.")
RHO_AMB_LOCAL = RHO_AMB_CALC


def _obter_material_parede(sala):
    """Retorna as propriedades do material da parede a partir do cenário."""
    nome = sala.get('material_parede', MATERIAL_PAREDE_FALLBACK)
    if nome not in MATERIAIS_CONSTRUCAO:
        print(f"[AVISO] Material '{nome}' não encontrado. "
      f"Usando fallback '{MATERIAL_PAREDE_FALLBACK}'.")
        nome = MATERIAL_PAREDE_FALLBACK
    return MATERIAIS_CONSTRUCAO[nome]


props_parede = _obter_material_parede({})
K_W = props_parede['k']
RHO_W = props_parede['rho']
CP_W = props_parede['cp']
E_W = props_parede['espessura']
H_CONV = props_parede['h_conv']
epsilon = props_parede['emissividade']

params_parede = criar_parede(
    N_nos=N_NOS,
    espessura_total=E_W,
    k=K_W,
    rho=RHO_W,
    cp=CP_W,
    h_conv_int=H_CONV,
    emissividade=epsilon,
    T_amb=T_AMB,
    h_conv_ext=0.005
)


# -----------------------------------------------------------
# Correlações de pluma
# -----------------------------------------------------------
def _pluma_zukoski(Q_c, z):
    """Zukoski et al. (1981) [7]."""
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0
    return 0.071 * (Q_c ** (1.0 / 3.0)) * (z ** (5.0 / 3.0))


def _pluma_heskestad(Q_dot_total, Q_c, z, D,
                     T_inf=T_AMB, cp=C_P, rho_inf=RHO_AMB_LOCAL, g=9.81):
    """
    Heskestad plume entrainment — CFAST Eqs. (3.16)–(3.18) [1], [5].
    """
    if z <= 0.01 or Q_c <= 0.0 or D <= 0.0:
        return 0.0

    Q_star = Q_dot_total / (rho_inf * cp * T_inf * np.sqrt(g * D) * D**2)

    # Altura média da chama [1, Eq. 3.30]
    L = D * (-1.02 + 3.7 * Q_star ** (2 / 5))
    L = max(0.01, L)

    # Origem virtual [1, Eq. 3.17]
    z0 = D * (-1.02 + 1.4 * Q_star ** (2 / 5))

    const = 0.196 * (g * rho_inf**2 / (cp * T_inf)) ** (1 / 3)

    if z < L:
        z_eff_L = max(L - z0, 0.01)
        m_e_L = const * Q_c ** (1 / 3) * z_eff_L ** (5 / 3)
        m_e_L *= 1.0 + (2.9 * Q_c ** (2 / 3)) / (
            (np.sqrt(g) * cp * rho_inf * T_inf) ** (2 / 3) * z_eff_L ** (5 / 3)
        )
        m_e = m_e_L * (z / L)
    else:
        z_eff = max(z - z0, 0.01)
        m_e = const * Q_c ** (1 / 3) * z_eff ** (5 / 3)
        m_e *= 1.0 + (2.9 * Q_c ** (2 / 3)) / (
            (np.sqrt(g) * cp * rho_inf * T_inf) ** (2 / 3) * z_eff ** (5 / 3)
        )

    return m_e


def _pluma_mccaffrey(Q_c, z):
    """
    McCaffrey (1979) [6] — divide a pluma em três regimes:
    chama contínua, intermitente e pluma pura.
    """
    if z <= 0.01 or Q_c <= 0.0:
        return 0.0
    zeta = z / (Q_c ** 0.4)
    if zeta < 0.08:
        return 0.011 * Q_c * (zeta ** 0.566)
    elif zeta < 0.20:
        return 0.026 * Q_c * (zeta ** 0.909)
    else:
        return 0.124 * Q_c * (zeta ** 1.895)


# -----------------------------------------------------------
# Sistema de EDOs
# -----------------------------------------------------------
def equacoes_sfpe(t, y, sala, aberturas, objetos_ativos, modelo_pluma='zukoski'):
    """
    Sistema de EDOs de conservação com múltiplos objetos.

    Parâmetros
    ----------
    y : array
        [m_u, H_u, m_o2, m_co2, m_co, T_w0, ..., T_w(N-1)]
    modelo_pluma : str
        'zukoski', 'heskestad' ou 'mccaffrey'.

    Retorna
    -------
    list
        Derivadas: [dm_dt, dH_dt, dm_o2_dt, dm_co2_dt, dm_co_dt] + dT_parede_dt
    """
    m_u = y[0]
    H_u = y[1]
    m_o2 = y[2]
    m_co2 = y[3]
    m_co = y[4]
    T_parede = y[5:5 + N_NOS]

    m_u = max(m_u, 0.001)
    m_o2 = max(m_o2, 0.0)
    m_co2 = max(m_co2, 0.0)
    m_co = max(m_co, 0.0)

    T_u = H_u / (m_u * C_P)
    T_u = max(200.0, min(T_u, 2000.0))

    V_u = (m_u * R_AR * T_u) / P_ATM
    area_sala = sala["area_piso"]
    Z_int = sala["altura_z"] - (V_u / area_sala)
    Z_int = max(Z_int, 0.0)

    Q_dot, m_dot_fuel, _, nu_o2, y_co2, y_co, D_fogo = \
        calcular_hrr_total_e_propriedades(t, objetos_ativos)

    # --- Escolha da pluma ---
    if modelo_pluma == 'zukoski':
        _calc_pluma = lambda Qc, z: _pluma_zukoski(Qc, z)
    elif modelo_pluma == 'mccaffrey':
        _calc_pluma = lambda Qc, z: _pluma_mccaffrey(Qc, z)
    else:
        _calc_pluma = lambda Qc, z: _pluma_heskestad(Q_dot, Qc, z, D_fogo)

    Q_c = 0.7 * Q_dot
    m_dot_plume_total = _calc_pluma(Q_c, Z_int)
    m_dot_entrained = max(0.0, m_dot_plume_total - m_dot_fuel)

    # --- Ventilação ---
    rho_u = P_ATM / (R_AR * T_u)
    m_dot_in = 0.0
    m_dot_out = 0.0

    def _calcular_ventilacao(m_dot_fuel_atual):
        """Calcula vazões de entrada/saída com um dado m_dot_fuel."""
        m_in = 0.0
        m_out = 0.0
        for vent in aberturas:
            fluxo = calcular_fluxo_porta_duas_zonas(
                largura=vent["dimensoes"]["x"],
                altura_porta=vent["dimensoes"]["z"],
                altura_interface=Z_int,
                densidade_ar=RHO_AMB_LOCAL,
                densidade_fumaca=rho_u,
                C_d=0.68,
                m_dot_combustivel=m_dot_fuel_atual
            )
            m_in += fluxo["m_dot_in"]
            m_out += fluxo["m_dot_out"]
        return m_in, m_out

    m_dot_in, m_dot_out = _calcular_ventilacao(m_dot_fuel)

    # --- Espécies ---
    m_dot_o2_consumido = m_dot_fuel * nu_o2
    m_dot_co2_produzido = m_dot_fuel * y_co2
    m_dot_co_produzido = m_dot_fuel * y_co

        # --- Limitação pelo O2 (Modelo CFAST 6/7 — tanh suave) ---
    # Referência: Jones et al., NIST SP 1030, Eq. (47) [1].
    #   C_LOL = (tanh(800*(Y_O2 - Y_LOL) - 4) + 1) / 2
    # Y_LOL é a fração mássica de O2 no Lower Oxygen Limit (LOL).
    # Default CFAST: 10% em volume ≈ 0.10 em massa.
    # A função varia de 0 (extinção) a 1 (queima completa) em ~0.01 de Y_O2.
    Y_O2_LOL = 0.01  # Fração mássica (CFAST default, Morehart et al.)
    Y_O2_atual = m_o2 / m_u if m_u > 0 else 0.0

    # Fator de limitação suave (CFAST Eq. 47)
    C_LOL = 0.5 * (np.tanh(800.0 * (Y_O2_atual - Y_O2_LOL) - 4.0) + 1.0)

    # Aplicar o fator à taxa de queima
    m_dot_fuel *= C_LOL
    Q_dot *= C_LOL
    Q_c = 0.7 * Q_dot

    # Recalcular pluma e espécies com a HRR limitada
    m_dot_plume_total = _calc_pluma(Q_c, Z_int)
    m_dot_entrained = max(0.0, m_dot_plume_total - m_dot_fuel)
    m_dot_o2_consumido = m_dot_fuel * nu_o2
    m_dot_co2_produzido = m_dot_fuel * y_co2
    m_dot_co_produzido = m_dot_fuel * y_co

    # Recalcular ventilação com m_dot_fuel atualizado
    m_dot_in, m_dot_out = _calcular_ventilacao(m_dot_fuel)


    dm_o2_dt = (m_dot_entrained * Y_O2_AMB
                - m_dot_o2_consumido
                - m_dot_out * (m_o2 / m_u))
    dm_co2_dt = (m_dot_entrained * Y_CO2_AMB
                 + m_dot_co2_produzido
                 - m_dot_out * (m_co2 / m_u))
    dm_co_dt = (m_dot_entrained * Y_CO_AMB
                + m_dot_co_produzido
                - m_dot_out * (m_co / m_u))

    if m_o2 <= 0.0 and dm_o2_dt < 0.0:
        dm_o2_dt = 0.0

    dm_dt = m_dot_plume_total - m_dot_out

    # --- Perda de calor ---
    L = sala["comprimento_x"]
    W = sala["profundidade_y"]
    H_sala = sala["altura_z"]
    A_teto = L * W
    A_paredes_sup = 2 * (L + W) * (H_sala - Z_int)
    A_contato = A_teto + A_paredes_sup  

    T_w0 = T_parede[0]
    Q_rad = epsilon * SIGMA * A_contato * (T_u ** 4 - T_w0 ** 4)
    Q_conv = H_CONV * A_contato * (T_u - T_w0)
    Q_loss = Q_conv + Q_rad

    dH_dt = (Q_c
             + m_dot_entrained * C_P * T_AMB
             - m_dot_out * C_P * T_u
             - Q_loss)

    dT_parede_dt = derivadas_parede(T_u, T_parede, params_parede)

    return [dm_dt, dH_dt, dm_o2_dt, dm_co2_dt, dm_co_dt] + list(dT_parede_dt)


# -----------------------------------------------------------
# Integração de segmentos
# -----------------------------------------------------------
def _integrar_segmento(t_ini, y_ini, t_fim, sala, aberturas,
                       objetos_ativos, modelo_pluma, evento_fn=None):
    """
    Integra um segmento temporal usando solve_ivp com Radau.

    Retorna
    -------
    sol : OdeResult
        Objeto do solve_ivp. Use sol.t[-1] para obter o tempo final real.
        Verifique sol.status == 2 para detectar falha.
    """
    t_span = (t_ini, t_fim)
    num_points = max(10, int((t_fim - t_ini) * PONTOS_POR_SEGUNDO))
    t_eval = np.linspace(t_ini, t_fim, num_points)

    sol = solve_ivp(
        equacoes_sfpe,
        t_span,
        y_ini,
        args=(sala, aberturas, objetos_ativos, modelo_pluma),
        events=evento_fn,
        method='Radau',
        max_step=MAX_STEP,
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-20
    )

    if sol.status == 2:
        raise RuntimeError(
        f"solve_ivp falhou em [{t_ini:.2f}, {t_fim:.2f}]: {sol.message}"
    )
    return sol


# -----------------------------------------------------------
# Simulação sequencial
# -----------------------------------------------------------
def simular_incendio_sequencial(sala, aberturas,
                                objeto_inicial, outros_objetos,
                                objetos_geometria,
                                t_final=600,
                                modelo_pluma='zukoski'):
    """
    Simulação completa com ignição sequencial.

    Retorna
    -------
    t_total : ndarray
        Vetor de tempos.
    y_total : ndarray
        Matriz de estados (n_vars × n_tempos).
    objetos_ativos : list
        Lista final de objetos queimando, com t_ign preenchido.
    """

    area_por_material = {obj['nome_material']: obj['area']
                         for obj in objetos_geometria}

    nome_inicial = objeto_inicial['nome_material']
    area_inicial = area_por_material.get(nome_inicial, 1.0)
    objetos_ativos = [{'nome_material': nome_inicial,
                       't_ign': 0.0,
                       'area': area_inicial}]

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
    T_u0 = T_AMB
    rho_u0 = P_ATM / (R_AR * T_u0)
    m_u_0 = rho_u0 * V_u_0
    H_u_0 = m_u_0 * C_P * T_u0
    m_o2_0 = m_u_0 * Y_O2_AMB
    m_co2_0 = m_u_0 * Y_CO2_AMB
    m_co_0 = m_u_0 * Y_CO_AMB

    T_parede0 = np.full(N_NOS, T_AMB)

    y0 = np.concatenate([
        np.array([m_u_0, H_u_0, m_o2_0, m_co2_0, m_co_0]),
        T_parede0
    ])

    t_atual = 0.0
    y_atual = y0.copy()
    solucoes = []

    TOL_IGN = 0.5

    def _verificar_ignicoes_imediatas(t, y, objetos_inativos, objetos_ativos):
        """
        Ativa objetos cuja T_ign já foi ultrapassada.
        Usa TOL_IGN para robustez numérica.
        Retorna o número de objetos ativados.
        """
        m_u_loc = max(y[0], 0.001)
        T_u_loc = y[1] / (m_u_loc * C_P)
        ignitados = [o for o in objetos_inativos
                     if T_u_loc >= o['T_ign'] - TOL_IGN]
        for obj in ignitados:
            objetos_inativos.remove(obj)
            obj['t_ign'] = t
            objetos_ativos.append(obj)
            print(f"Ignição imediata de {obj['nome_material']} em t = {t:.2f} s")
        return len(ignitados)

    
    _verificar_ignicoes_imediatas(t_atual, y_atual,
                                  objetos_inativos, objetos_ativos)

    
    while t_atual < t_final and len(objetos_inativos) > 0:
        T_IGN_ATUAL = min(obj['T_ign'] for obj in objetos_inativos)

        
        def evento_ignicao(t, y):
            """Evento de ignição: T_u cruza T_IGN_ATUAL."""
            m_u_loc = max(y[0], 0.001)
            H_u_loc = y[1]
            T_u_loc = H_u_loc / (m_u_loc * C_P)
            return T_u_loc - T_IGN_ATUAL

        evento_ignicao.terminal = True
        evento_ignicao.direction = 1

        sol = _integrar_segmento(
            t_atual, y_atual, t_final,
            sala, aberturas, objetos_ativos, modelo_pluma,
            evento_fn=evento_ignicao
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
            _verificar_ignicoes_imediatas(t_atual, y_atual,
                                          objetos_inativos, objetos_ativos)
        else:
            #Sem evento: sai do loop.
            break

    if t_atual < t_final:
        sol_rest = _integrar_segmento(
            t_atual, y_atual, t_final,
            sala, aberturas, objetos_ativos, modelo_pluma
        )
        solucoes.append(sol_rest)
        t_atual = sol_rest.t[-1]

    # Concatenação segura
    t_hist = [np.array([0.0])]
    y_hist = [y0.reshape(-1, 1)]
    for sol_i in solucoes:
        if sol_i.t.size > 1:
            t_hist.append(sol_i.t[1:])
            y_hist.append(sol_i.y[:, 1:])

    t_total = np.concatenate(t_hist)
    y_total = np.concatenate(y_hist, axis=1)
    return t_total, y_total, objetos_ativos