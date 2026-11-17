"""
Módulo de Física do Simulador de Incêndio
Responsável por calcular fluxos de ventilação 
e resolver o balanço de massa (Fluidodinâmica).

Modelo de escoamento por aberturas verticais (Emmons / Rockett), com a
camada inferior à densidade ambiente. Todas as aberturas de um mesmo
compartimento compartilham UMA única diferença de pressão de referência,
o que equivale a um único plano neutro Z_N (altura absoluta, medida do
piso), resolvido de modo que  Σ m_out − Σ m_in = m_dot_combustivel.
"""

import numpy as np
from scipy.optimize import brentq

G = 9.81


def _fluxos_abertura(largura: float, z_base: float, z_topo: float,
                     Z_i: float, Z_N: float,
                     densidade_ar: float, densidade_fumaca: float,
                     C_d: float):
    """
    Vazões através de UMA abertura vertical, dado o plano neutro Z_N
    (altura absoluta, com Z_N >= Z_i).

    Perfil de pressão interior − exterior (drho = rho_ar − rho_fumaça):
        z < Z_i         : Δp = −drho·g·(Z_N − Z_i)   (constante → entrada)
        Z_i <= z < Z_N  : Δp = −drho·g·(Z_N − z)     (entrada)
        z >= Z_N        : Δp = +drho·g·(z − Z_N)     (saída de fumaça)

    Retorna (m_in_inf, m_in_sup, m_out) [kg/s]:
        m_in_inf : entrada abaixo da interface
        m_in_sup : entrada entre a interface e o plano neutro
        m_out    : saída acima do plano neutro
    """
    drho = densidade_ar - densidade_fumaca
    C_in  = C_d * largura * np.sqrt(2.0 * densidade_ar    * G * drho)
    C_out = C_d * largura * np.sqrt(2.0 * densidade_fumaca * G * drho)

    # Entrada abaixo da interface (Δp constante)
    h_inf = max(0.0, min(z_topo, Z_i) - z_base)
    m_in_inf = C_in * np.sqrt(Z_N - Z_i) * h_inf

    # Entrada entre a interface e o plano neutro (Δp linear)
    a = max(z_base, Z_i)
    b = min(z_topo, Z_N)
    m_in_sup = 0.0
    if b > a:
        m_in_sup = C_in * (2.0 / 3.0) * ((Z_N - a) ** 1.5 - (Z_N - b) ** 1.5)

    # Saída acima do plano neutro (Δp linear)
    a_out = max(z_base, Z_N)
    m_out = 0.0
    if z_topo > a_out:
        m_out = C_out * (2.0 / 3.0) * ((z_topo - Z_N) ** 1.5
                                       - (a_out - Z_N) ** 1.5)

    return m_in_inf, m_in_sup, m_out


def calcular_fluxos_compartimento(
    aberturas: list,
    altura_interface: float,
    densidade_ar: float,
    densidade_fumaca: float,
    C_d: float = 0.68,
    m_dot_combustivel: float = 0.0,
    m_dot_expansao: float = 0.0,
    razao_T_out: float = 1.0,
) -> dict:
    """
    Resolve o escoamento de TODAS as aberturas de um compartimento com um
    único plano neutro Z_N impondo o fechamento global quase-isobárico
    (conservação volumétrica dV_total/dt = 0 / Lei dos Gases Ideais).

    aberturas      : lista de dicts com 'largura', 'altura' e 'z_base' [m]
    m_dot_expansao : termo líquido equivalente de geração térmica
                     (Q_c - Q_loss) / (c_p * T_amb) [kg/s]
    razao_T_out    : razão T_u / T_amb (= rho_ar / rho_u) para ponderar a
                     expansão específica na exaustão; com 1.0 recupera
                     Σm_out − Σm_in − m_fuel − m_dot_expansao = 0

    Retorna dict com:
        Z_N, m_dot_in, m_dot_in_inf, m_dot_in_sup, m_dot_out, m_dot_exp (totais) e
        'por_abertura': lista de dicts (z_n relativo à base, m_dot_in,
        m_dot_in_inf, m_dot_in_sup, m_dot_out) na mesma ordem de entrada.
    """
    Z_i = altura_interface
    geo = [(ab["largura"], ab["z_base"], ab["z_base"] + ab["altura"])
           for ab in aberturas
           if ab["largura"] > 0.0 and ab["altura"] > 0.0]

    def _montar(Z_N, fluxos_ab):
        por_abertura = []
        tot_inf, tot_sup, tot_out = 0.0, 0.0, 0.0
        for (W, zb, zt), (m_inf, m_sup, m_out) in zip(geo, fluxos_ab):
            por_abertura.append({
                "z_n": max(0.0, min(zt - zb, Z_N - zb)),
                "m_dot_in": m_inf + m_sup,
                "m_dot_in_inf": m_inf,
                "m_dot_in_sup": m_sup,
                "m_dot_out": m_out,
            })
            tot_inf += m_inf
            tot_sup += m_sup
            tot_out += m_out
        m_exp_liquido = m_dot_expansao - (razao_T_out - 1.0) * tot_out
        return {"Z_N": Z_N,
                "m_dot_in": tot_inf + tot_sup,
                "m_dot_in_inf": tot_inf,
                "m_dot_in_sup": tot_sup,
                "m_dot_out": tot_out,
                "m_dot_exp": m_exp_liquido,
                "por_abertura": por_abertura}

    # Sem estratificação ou sem aberturas válidas: não há escoamento
    Z_topo_max = max((zt for _, _, zt in geo), default=0.0)
    if densidade_fumaca >= densidade_ar or not geo or Z_i >= Z_topo_max:
        return _montar(Z_i, [(0.0, 0.0, 0.0)] * len(geo))

    def _avaliar(Z_N):
        return [_fluxos_abertura(W, zb, zt, Z_i, Z_N,
                                 densidade_ar, densidade_fumaca, C_d)
                for W, zb, zt in geo]

    def f(Z_N: float) -> float:
        fl = _avaliar(Z_N)
        m_out_total = sum(m_out for _, _, m_out in fl)
        m_in_total = sum(m_inf + m_sup for m_inf, m_sup, _ in fl)
        return (razao_T_out * m_out_total
                - m_in_total
                - m_dot_combustivel
                - m_dot_expansao)

    # f é estritamente decrescente em [Z_i, Z_topo_max]:
    # Se f(Z_i) <= 0 (engasgamento / choking), plano neutro na interface.
    # Se f(Z_topo_max) >= 0 (contração térmica no resfriamento onde m_dot_exp < 0),
    # plano neutro no topo da abertura mais alta mantendo confinamento seguro para brentq.
    f_i = f(Z_i)
    if f_i <= 0.0:
        return _montar(Z_i, _avaliar(Z_i))

    f_topo = f(Z_topo_max)
    if f_topo >= 0.0:
        return _montar(Z_topo_max, _avaliar(Z_topo_max))

    Z_N = brentq(f, Z_i, Z_topo_max, xtol=1e-12, rtol=1e-12, maxiter=200)
    return _montar(Z_N, _avaliar(Z_N))


def calcular_fluxo_porta_duas_zonas(
    largura: float,
    altura_porta: float,
    altura_interface: float,
    densidade_ar: float,
    densidade_fumaca: float,
    C_d: float = 0.68,
    m_dot_combustivel: float = 0.0,
    z_base_porta: float = 0.0,
    m_dot_expansao: float = 0.0,
    razao_T_out: float = 1.0,
) -> dict:
    """
    Interface de compatibilidade para UMA abertura isolada.

    ATENÇÃO: para compartimentos com mais de uma abertura use
    calcular_fluxos_compartimento(); chamar esta função por abertura
    impõe m_dot_combustivel e m_dot_expansao em cada uma (contagem múltipla).

    z_n é retornado relativo à base da abertura.
    """
    if densidade_fumaca >= densidade_ar:
        return {"z_n": altura_porta / 2.0, "m_dot_in": 0.0, "m_dot_out": 0.0,
                "m_dot_in_inf": 0.0, "m_dot_in_sup": 0.0}

    res = calcular_fluxos_compartimento(
        [{"largura": largura, "altura": altura_porta, "z_base": z_base_porta}],
        altura_interface, densidade_ar, densidade_fumaca,
        C_d=C_d, m_dot_combustivel=m_dot_combustivel,
        m_dot_expansao=m_dot_expansao, razao_T_out=razao_T_out,
    )
    if not res["por_abertura"]:
        return {"z_n": altura_porta, "m_dot_in": 0.0, "m_dot_out": 0.0,
                "m_dot_in_inf": 0.0, "m_dot_in_sup": 0.0}
    return dict(res["por_abertura"][0])
