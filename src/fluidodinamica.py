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
    z_base_porta: float = 0.0,
) -> dict:
    if densidade_fumaca >= densidade_ar:
        return {"z_n": altura_porta / 2.0, "m_dot_in": 0.0, "m_dot_out": 0.0}

    g = 9.81
    W = largura
    H = altura_porta
    drho = densidade_ar - densidade_fumaca

    z_i_raw = altura_interface - z_base_porta

    if z_i_raw >= H:
        return {"z_n": H, "m_dot_in": 0.0, "m_dot_out": 0.0}

    z_i = max(0.0, z_i_raw)

    C_in  = C_d * W * np.sqrt(2.0 * densidade_ar    * g * drho)
    C_out = (2.0 / 3.0) * C_d * W * np.sqrt(2.0 * densidade_fumaca * g * drho)

    def fluxos(z_n: float):
        if z_n < z_i or z_n > H:
            return np.inf, np.inf
        dz = z_n - z_i
        m_in  = C_in  * (z_i * np.sqrt(dz) + (2.0 / 3.0) * dz ** 1.5)
        m_out = C_out * (H - z_n) ** 1.5
        return m_in, m_out

    m_out_max = C_out * (H - z_i) ** 1.5
    if m_dot_combustivel >= m_out_max:
        return {"z_n": z_i, "m_dot_in": 0.0, "m_dot_out": m_out_max}

    def f(z: float) -> float:
        m_in, m_out = fluxos(z)
        if not (np.isfinite(m_in) and np.isfinite(m_out)):
            return np.nan
        return m_out - m_in - m_dot_combustivel

    r = densidade_fumaca / densidade_ar
    z0 = H * r ** (1.0 / 3.0) / (1.0 + r ** (1.0 / 3.0))
    z0 = max(z_i + 1e-6, min(z0, H - 1e-6))

    try:
        z_n_sol = newton(f, z0, tol=1e-8, maxiter=50)
        z_n_sol = max(z_i, min(H, z_n_sol))
    except (RuntimeError, ValueError):
        z_n_sol = H / (1.0 + (densidade_ar / densidade_fumaca) ** (1.0 / 3.0))
        z_n_sol = max(z_i, min(H, z_n_sol))

    m_in, m_out = fluxos(z_n_sol)
    return {"z_n": z_n_sol, "m_dot_in": m_in, "m_dot_out": m_out}