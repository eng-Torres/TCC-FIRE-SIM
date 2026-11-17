# src/edos_solver.py
"""
Módulo de Dinâmica de Camadas
Resolve as Equações de Conservação de Massa, Energia e Espécies.
Baseado no SFPE Handbook (Cap. 36) e CFAST (NIST SP 1026/1030).
Inclui correlação da pluma (Zukoski, McCaffrey ou Heskestad),
modelo de parede com condução transiente (N nós),
suporte a múltiplos objetos com ignição sequencial.
"""

import numpy as np
from scipy.integrate import solve_ivp

from src.fluidodinamica import calcular_fluxos_compartimento
from src.combustao import calcular_hrr_total_e_propriedades
from src.constantes import *
from src.banco_construcao import MATERIAIS_CONSTRUCAO
from src.termica_parede import criar_parede, derivadas_parede, calcular_perda_calor_superior

# Material de parede usado quando o cenário não especifica
MATERIAL_PAREDE_FALLBACK = 'gesso_16mm'

# Sanity check de RHO_AMB vs P_ATM/(R_AR*T_AMB)
RHO_AMB_CALC = P_ATM / (R_AR * T_AMB)
if abs(RHO_AMB - RHO_AMB_CALC) / RHO_AMB_CALC > 1e-3:
    print(f"[AVISO] RHO_AMB ({RHO_AMB:.4f}) difere de "
          f"P_ATM/(R_AR*T_AMB) ({RHO_AMB_CALC:.4f}). Usando valor calculado.")
RHO_AMB_LOCAL = RHO_AMB_CALC


def _obter_material_parede(sala):
    """Retorna as propriedades do material da parede a partir do cenário."""
    nome = sala.get('material_parede', MATERIAL_PAREDE_FALLBACK)
    if nome not in MATERIAIS_CONSTRUCAO:
        print(f"[AVISO] Material '{nome}' não encontrado. "
              f"Usando fallback '{MATERIAL_PAREDE_FALLBACK}'.")
        nome = MATERIAL_PAREDE_FALLBACK
    return MATERIAIS_CONSTRUCAO[nome]


def _construir_params_parede(sala, n_nos):
    """Cria o dicionário de parâmetros da parede para este cenário."""
    props = _obter_material_parede(sala)
    return criar_parede(
        N_nos=n_nos,
        espessura_total=props['espessura'],
        k=props['k'],
        rho=props['rho'],
        cp=props['cp'],
        h_conv_int=props['h_conv'],
        emissividade=props['emissividade'],
        T_amb=T_AMB,
        h_conv_ext=0.005,
    )


# -----------------------------------------------------------
# Correlações de pluma
# -----------------------------------------------------------
def _pluma_zukoski(Q_c, z):
    if z <= 0.0 or Q_c <= 0.0:
        return 0.0
    return 0.071 * (Q_c ** (1.0 / 3.0)) * (z ** (5.0 / 3.0))


def _pluma_heskestad(Q_dot_total, Q_c, z, D,
                     T_inf=T_AMB, cp=C_P, rho_inf=RHO_AMB_LOCAL, g=9.81):
    """Heskestad plume entrainment — CFAST Eqs. (3.16)-(3.18)."""
    if z <= 0.0 or Q_c <= 0.0 or D <= 0.0:
        return 0.0

    Q_star = Q_dot_total / (rho_inf * cp * T_inf * np.sqrt(g * D) * D**2)
    L = max(0.01, D * (-1.02 + 3.7 * Q_star ** (2 / 5)))
    # Origem virtual de Heskestad (1984): pode ser negativa (abaixo da base do fogo)
    # para fontes extensas / baixo Q* (Q* < 0.45), aumentando a distância efetiva (z - z0).
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
    if z <= 0.0 or Q_c <= 0.0:
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
def equacoes_sfpe(t, y, sala, aberturas, objetos_ativos, modelo_pluma,
                  n_nos, params_parede):
    """
    y = [m_u, H_u, m_o2, m_co2, m_co, T_w0, ..., T_w(N-1)]
    """
    m_u = max(y[0], 0.001)
    H_u = y[1]
    m_o2 = max(y[2], 0.0)
    m_co2 = max(y[3], 0.0)
    m_co = max(y[4], 0.0)
    T_parede = y[5:5 + n_nos]

    T_u = max(200.0, min(H_u / (m_u * C_P), 2000.0))

    V_u = (m_u * R_AR * T_u) / P_ATM
    area_sala = sala["area_piso"]
    Z_int = max(sala["altura_z"] - (V_u / area_sala), 0.0)

    Q_dot, m_dot_fuel, _, nu_o2, y_co2, y_co, D_fogo = \
        calcular_hrr_total_e_propriedades(t, objetos_ativos)

    if modelo_pluma == 'zukoski':
        _calc_pluma = lambda Qc, z: _pluma_zukoski(Qc, z)
    elif modelo_pluma == 'mccaffrey':
        _calc_pluma = lambda Qc, z: _pluma_mccaffrey(Qc, z)
    else:
        _calc_pluma = lambda Qc, z: _pluma_heskestad(Q_dot, Qc, z, D_fogo)

    rho_u = P_ATM / (R_AR * T_u)

    def _calcular_ventilacao(m_dot_fuel_atual, m_dot_exp_atual):
        """
        Retorna (m_in, m_in_sup, m_out):
            m_in     : entrada total de ar ambiente [kg/s]
            m_in_sup : parcela da entrada que ocorre ACIMA da interface [kg/s]
            m_out    : saída de fumaça da camada superior [kg/s]
        """
        # Todas as aberturas do compartimento compartilham UM plano neutro;
        # o combustível e a expansão térmica entram uma única vez no fechamento.
        geo_aberturas = [{
            "largura": vent["dimensoes"]["x"],
            "altura": vent["dimensoes"]["z"],
            "z_base": vent.get("posicao", {}).get("z", 0.0),
        } for vent in aberturas]
        fluxo = calcular_fluxos_compartimento(
            geo_aberturas,
            altura_interface=Z_int,
            densidade_ar=RHO_AMB_LOCAL,
            densidade_fumaca=rho_u,
            C_d=0.68,
            m_dot_combustivel=m_dot_fuel_atual,
            m_dot_expansao=m_dot_exp_atual,
            razao_T_out=T_u / T_AMB,
        )
        return fluxo["m_dot_in"], fluxo["m_dot_in_sup"], fluxo["m_dot_out"]

    # Limitação por O2 (CFAST 6/7, tanh suave)
    Y_O2_LOL = 0.01
    Y_O2_atual = m_o2 / m_u if m_u > 0 else 0.0
    C_LOL = 0.5 * (np.tanh(800.0 * (Y_O2_atual - Y_O2_LOL) - 4.0) + 1.0)

    m_dot_fuel *= C_LOL
    Q_dot *= C_LOL
    Q_c = 0.7 * Q_dot
    Q_r = 0.3 * Q_dot

    m_dot_plume_total = _calc_pluma(Q_c, Z_int)
    # Arrastamento de ar ambiente pela pluma (nunca negativo). O combustível
    # entra na camada SEMPRE, mesmo quando a pluma é curta (plume < fuel)
    # ou quando a interface atinge o piso (pluma = 0).
    m_dot_entrained = max(0.0, m_dot_plume_total - m_dot_fuel)
    m_dot_o2_consumido = m_dot_fuel * nu_o2
    m_dot_co2_produzido = m_dot_fuel * y_co2
    m_dot_co_produzido = m_dot_fuel * y_co

    # Perda de calor (pré-calculada para acoplar expansão térmica no plano neutro)
    T_w0 = T_parede[0]
    Q_loss, A_contato = calcular_perda_calor_superior(
        sala,
        Z_int,
        T_u,
        T_w0,
        params_parede['h_conv_int'],
        params_parede['emissividade'],
        sigma=SIGMA,
    )
    m_dot_exp_termica = (Q_c - Q_loss) / (C_P * T_AMB)

    m_dot_in, m_dot_in_sup, m_dot_out = _calcular_ventilacao(
        m_dot_fuel, m_dot_exp_termica
    )

    # Ar ambiente que chega à camada superior: apenas o arrastado pela pluma.
    # HIPÓTESE: o ar que entra pelas aberturas (m_dot_in, inclusive a parcela
    # m_dot_in_sup que entra acima da interface) é mais denso que a fumaça,
    # afunda e é depositado na camada inferior (como no CFAST). Como a camada
    # inferior é implícita (reservatório a T_AMB, Y_O2_AMB), essa massa não
    # entra nos balanços abaixo. Depositar m_dot_in_sup na camada superior
    # cria uma realimentação positiva espúria (interface desce -> mais ar
    # entra acima dela -> interface desce mais) que leva Z_int ao piso.
    # O tratamento completo exige uma camada inferior explícita.
    m_dot_ar_sup = m_dot_entrained

    dm_o2_dt = (m_dot_ar_sup * Y_O2_AMB
                - m_dot_o2_consumido
                - m_dot_out * (m_o2 / m_u))
    dm_co2_dt = (m_dot_ar_sup * Y_CO2_AMB
                 + m_dot_co2_produzido
                 - m_dot_out * (m_co2 / m_u))
    dm_co_dt = (m_dot_ar_sup * Y_CO_AMB
                + m_dot_co_produzido
                - m_dot_out * (m_co / m_u))

    if m_o2 <= 0.0 and dm_o2_dt < 0.0:
        dm_o2_dt = 0.0

    # Balanço de massa com as MESMAS parcelas usadas nas espécies e na energia
    dm_dt = m_dot_ar_sup + m_dot_fuel - m_dot_out

    # Toda massa que entra na camada traz sua entalpia (ar e combustível
    # a T_AMB); a fumaça sai a T_u.
    dH_dt = (Q_c
             + (m_dot_ar_sup + m_dot_fuel) * C_P * T_AMB
             - m_dot_out * C_P * T_u
             - Q_loss)

    # A parcela radiativa da chama (Q_r = 0.3 * Q_dot) incide como fluxo
    # adicional na superfície interna da parede (nó 0), conservando a energia
    # total liberada na combustão (Q_c + Q_r = Q_dot).
    q_rad_chama = Q_r / A_contato
    dT_parede_dt = derivadas_parede(T_u, T_parede, params_parede)
    dT_parede_dt[0] += q_rad_chama / params_parede['cap_surf']

    return [dm_dt, dH_dt, dm_o2_dt, dm_co2_dt, dm_co_dt] + list(dT_parede_dt)


# -----------------------------------------------------------
# Integração de segmentos
# -----------------------------------------------------------
def _integrar_segmento(t_ini, y_ini, t_fim, sala, aberturas,
                       objetos_ativos, modelo_pluma,
                       n_nos, params_parede,
                       max_step, pontos_por_segundo,
                       evento_fn=None):
    t_span = (t_ini, t_fim)
    num_points = max(10, int((t_fim - t_ini) * pontos_por_segundo))
    t_eval = np.linspace(t_ini, t_fim, num_points)

    sol = solve_ivp(
        equacoes_sfpe,
        t_span,
        y_ini,
        args=(sala, aberturas, objetos_ativos, modelo_pluma,
              n_nos, params_parede),
        events=evento_fn,
        method='Radau',
        max_step=max_step,
        t_eval=t_eval,
        rtol=1e-6,
        atol=1e-9,
    )
    if sol.status < 0:
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
                                modelo_pluma='zukoski',
                                max_step=0.1,
                                pontos_por_segundo=10,
                                n_nos=5):
    """
    Simulação completa com ignição sequencial.

    Retorna
    -------
    t_total, y_total, objetos_ativos
    """
    if n_nos < 2:
        raise ValueError(f"n_nos deve ser >= 2. Recebido: {n_nos}")

    # Reconstrói params_parede para este cenário (material + n_nos)
    params_parede_local = _construir_params_parede(sala, n_nos)

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
            'area': area,
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

    T_parede0 = np.full(n_nos, T_AMB)
    y0 = np.concatenate([
        np.array([m_u_0, H_u_0, m_o2_0, m_co2_0, m_co_0]),
        T_parede0,
    ])

    t_atual = 0.0
    y_atual = y0.copy()
    solucoes = []

    TOL_IGN = 0.5

    def _verificar_ignicoes_imediatas(t, y, inativos, ativos):
        m_u_loc = max(y[0], 0.001)
        T_u_loc = y[1] / (m_u_loc * C_P)
        ignitados = [o for o in inativos
                     if T_u_loc >= o['T_ign'] - TOL_IGN]
        for obj in ignitados:
            inativos.remove(obj)
            obj['t_ign'] = t
            ativos.append(obj)
            print(f"Ignição imediata de {obj['nome_material']} em t = {t:.2f} s")
        return len(ignitados)

    _verificar_ignicoes_imediatas(t_atual, y_atual,
                                  objetos_inativos, objetos_ativos)

    while t_atual < t_final and len(objetos_inativos) > 0:
        T_IGN_ATUAL = min(obj['T_ign'] for obj in objetos_inativos)

        def evento_ignicao(t, y, *args):
            m_u_loc = max(y[0], 0.001)
            T_u_loc = y[1] / (m_u_loc * C_P)
            return T_u_loc - T_IGN_ATUAL

        evento_ignicao.terminal = True
        evento_ignicao.direction = 1

        sol = _integrar_segmento(
            t_atual, y_atual, t_final,
            sala, aberturas, objetos_ativos, modelo_pluma,
            n_nos, params_parede_local,
            max_step, pontos_por_segundo,
            evento_fn=evento_ignicao,
        )
        solucoes.append(sol)

        if sol.t_events[0].size > 0 and len(objetos_inativos) > 0:
            t_atual = float(sol.t_events[0][0])
            y_atual = sol.y_events[0][0].copy()
            if sol.t.size == 0 or sol.t[-1] < t_atual:
                sol.t = np.append(sol.t, t_atual)
                sol.y = np.column_stack([sol.y, y_atual])

            obj_ignitado = min(objetos_inativos, key=lambda x: x['T_ign'])
            objetos_inativos.remove(obj_ignitado)
            obj_ignitado['t_ign'] = t_atual
            objetos_ativos.append(obj_ignitado)
            print(f"Ignição de {obj_ignitado['nome_material']} em t = {t_atual:.2f} s")
            _verificar_ignicoes_imediatas(t_atual, y_atual,
                                          objetos_inativos, objetos_ativos)
        else:
            t_atual = sol.t[-1]
            y_atual = sol.y[:, -1]
            break

    if t_atual < t_final:
        sol_rest = _integrar_segmento(
            t_atual, y_atual, t_final,
            sala, aberturas, objetos_ativos, modelo_pluma,
            n_nos, params_parede_local,
            max_step, pontos_por_segundo,
        )
        solucoes.append(sol_rest)
        t_atual = sol_rest.t[-1]

    t_hist = [np.array([0.0])]
    y_hist = [y0.reshape(-1, 1)]
    for sol_i in solucoes:
        if sol_i.t.size > 1:
            t_hist.append(sol_i.t[1:])
            y_hist.append(sol_i.y[:, 1:])

    t_total = np.concatenate(t_hist)
    y_total = np.concatenate(y_hist, axis=1)
    return t_total, y_total, objetos_ativos
