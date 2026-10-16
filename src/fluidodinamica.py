"""
Módulo de Física do Simulador de Incêndio
Responsável por calcular fluxos de ventilação 
e resolver o balanço de massa (Fluidodinâmica).
"""

import numpy as np
from scipy.optimize import newton

def calcular_fluxo_porta_duas_zonas(
    largura: float,
    altura_porta: float,
    altura_interface: float,
    densidade_ar: float,
    densidade_fumaca: float,
    C_d: float = 0.68,
    m_dot_combustivel: float = 0.0,
    z_base_porta: float=0.0,
) -> dict:
    if densidade_fumaca >= densidade_ar:
        return {"z_n": altura_porta / 2.0, "m_dot_in": 0.0, "m_dot_out": 0.0}

    g = 9.81
    W = largura
    H = altura_porta
    z_i = max(0.0, min(H, altura_interface-z_base_porta))
    drho = densidade_ar - densidade_fumaca

    # Constantes auxiliares
    C_in  = C_d * W * np.sqrt(2.0 * densidade_ar * g * drho)
    C_out = (2.0/3.0) * C_d * W * np.sqrt(2.0 * densidade_fumaca * g * drho)

    def fluxos(z_n):
        if z_n < z_i or z_n > H:
            return np.inf, np.inf
        # Entrada: 0 → z_i (velocidade uniforme) + z_i → z_n (perfil triangular)
        termo1 = z_i * np.sqrt(z_n - z_i)
        termo2 = (2.0/3.0) * (z_n - z_i)**1.5
        m_in = C_in * (termo1 + termo2)
        # Saída: z_n → H (perfil triangular)
        m_out = C_out * (H - z_n)**1.5
        return m_in, m_out

    def f(z):
        m_in, m_out = fluxos(z)
        if np.isinf(m_in) or np.isinf(m_out):
            return np.inf
        return m_out - m_in - m_dot_combustivel

    # Chute inicial
    r = densidade_fumaca / densidade_ar
    z0 = H * r**(1/3) / (1 + r**(1/3))
    z0 = max(z_i + 1e-3, min(z0, H - 1e-3))

    try:
        z_n_sol = newton(f, z0, tol=1e-8, maxiter=50)
        z_n_sol = max(z_i, min(H, z_n_sol))
    except RuntimeError:
        z_n_sol = H / (1 + (densidade_ar/densidade_fumaca)**(1/3))
        z_n_sol = max(z_i, min(H, z_n_sol))

    m_in, m_out = fluxos(z_n_sol)
    return {"z_n": z_n_sol, "m_dot_in": m_in, "m_dot_out": m_out}