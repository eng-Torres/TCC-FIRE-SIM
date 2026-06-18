"""
Módulo de transferência de calor na parede (condução transiente 1D).
Discretização por diferenças finitas com N nós.
Baseado em Jones et al. (2009) e Incropera et al. (2007).
"""
import numpy as np
from src.constantes import SIGMA

def criar_parede(N_nos, espessura_total, k, rho, cp,
                 h_conv_int, emissividade, T_amb,
                 h_conv_ext=0.005):
    """
    Inicializa as constantes geométricas e térmicas da parede.

    Retorna um dicionário `params` que será usado pela função
    `derivadas_parede` em cada avaliação das EDOs.
    """
    dx = espessura_total / (N_nos - 1)          # distância entre nós [m]
    cap_node = rho * cp * dx                    # capacitância interior [kJ/(m²·K)]
    cap_surf = cap_node / 2.0                   # capacitância da superfície
    cond = k / dx                               # condutância entre nós [kW/(m²·K)]

    params = {
        'N': N_nos,
        'dx': dx,
        'cap_node': cap_node,
        'cap_surf': cap_surf,
        'cond': cond,
        'h_conv_int': h_conv_int,
        'emissividade': emissividade,
        'T_amb': T_amb,
        'h_conv_ext': h_conv_ext,
        'sigma': SIGMA
    }
    return params


def derivadas_parede(T_u, T_parede, params):
    """
    Calcula as derivadas temporais das temperaturas dos nós da parede.

    Argumentos:
        T_u        : temperatura da camada superior [K]
        T_parede   : array 1D com as temperaturas dos nós [K]
        params     : dicionário retornado por criar_parede()

    Retorna:
        dT_dt      : array 1D com as derivadas dT/dt para cada nó [K/s]
    """
    N = params['N']
    cond = params['cond']
    h_int = params['h_conv_int']
    eps = params['emissividade']
    sigma = params['sigma']
    T_amb = params['T_amb']
    h_ext = params['h_conv_ext']
    cap_surf = params['cap_surf']
    cap_node = params['cap_node']

    dT_dt = np.zeros(N)

    # ---- Nó da superfície interna (índice 0) ----
    q_in = h_int * (T_u - T_parede[0]) + eps * sigma * (T_u**4 - T_parede[0]**4)
    q_cond = cond * (T_parede[0] - T_parede[1])
    dT_dt[0] = (q_in - q_cond) / cap_surf

    # ---- Nós interiores ----
    for i in range(1, N-1):
        q_cond_entra = cond * (T_parede[i-1] - T_parede[i])
        q_cond_sai   = cond * (T_parede[i] - T_parede[i+1])
        dT_dt[i] = (q_cond_entra - q_cond_sai) / cap_node

    # ---- Nó da face externa (índice N-1) ----
    q_cond = cond * (T_parede[N-2] - T_parede[N-1])
    q_out  = h_ext * (T_parede[N-1] - T_amb)
    dT_dt[N-1] = (q_cond - q_out) / cap_surf

    return dT_dt