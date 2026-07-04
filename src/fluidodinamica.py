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
) -> dict:
    """
    Fluxos de massa bidirecionais em uma porta vertical – formulação correta.
    Hipótese: camada inferior = ar exterior (densidade_ar).
    """

    # Não há empuxo se a densidade da fumaça ≥ densidade do ar
    if densidade_fumaca >= densidade_ar:
        return {"z_n": altura_porta / 2.0, "m_dot_in": 0.0, "m_dot_out": 0.0}

    g = 9.81
    W = largura
    H = altura_porta
    z_i = max(0.0, min(altura_porta, altura_interface))
    drho = densidade_ar - densidade_fumaca
    rho_a = densidade_ar
    rho_f = densidade_fumaca

    # Diferença de pressão constante abaixo da interface (ΔP0)
    delta_P0 = g * drho * (H - z_i)   # positivo se T_u > T_amb

    # Se ΔP0 for negativo ou nulo, não há escoamento
    if delta_P0 <= 0.0:
        return {"z_n": H, "m_dot_in": 0.0, "m_dot_out": 0.0}

    # ------------------------------------------------------------
    # Apenas o Caso B é possível: plano neutro na camada superior (z_n ≥ z_i)
    # ------------------------------------------------------------
    def fluxos(z_n):
        """Calcula m_in e m_out para um dado z_n (z_n entre z_i e H)."""
        if z_n < z_i or z_n > H:
            return np.inf, np.inf

        # Entrada de ar frio (0 ≤ z ≤ z_i): velocidade uniforme
        m_in = C_d * W * z_i * np.sqrt(2.0 * rho_a * delta_P0)

        # Saída de fumaça (z_n ≤ z ≤ H): integração com perfil linear
        # ΔP(z) = delta_P0 + g*drho*(z - z_i)  → em z = z_n, ΔP=0
        # A integral ∫ √(g*drho*(z - z_n)) dz de z_n a H dá (2/3)*(H - z_n)^(3/2)
        m_out = (2.0 / 3.0) * C_d * W * np.sqrt(2.0 * rho_f * g * drho) * (H - z_n) ** 1.5

        return m_in, m_out

    # ------------------------------------------------------------
    # Função objetivo: m_in - m_out - m_dot_combustivel = 0
    # ------------------------------------------------------------
    def f(z):
        m_in, m_out = fluxos(z)
        if np.isinf(m_in) or np.isinf(m_out):
            return np.inf
        return m_in - m_out - m_dot_combustivel

    # Chute inicial para z_n (entre z_i e H)
    # Solução aproximada para camada homogênea (z_i = 0) como referência
    r = densidade_fumaca / densidade_ar
    z0 = H * r**(1/3) / (1 + r**(1/3))
    # Garantir que o chute está dentro do domínio [z_i, H]
    z0 = max(z_i + 0.01, min(z0, H - 0.01))

    try:
        z_n_sol = newton(f, z0, tol=1e-8, maxiter=50)
        z_n_sol = max(z_i, min(H, z_n_sol))
    except RuntimeError:
        # Fallback analítico (ignorando interface)
        z_n_sol = H / (1 + (densidade_ar / densidade_fumaca)**(1.0/3.0))
        z_n_sol = max(z_i, min(H, z_n_sol))

    m_in_final, m_out_final = fluxos(z_n_sol)
    return {"z_n": z_n_sol, "m_dot_in": m_in_final, "m_dot_out": m_out_final}