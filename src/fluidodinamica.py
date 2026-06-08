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
    Fluxos de massa bidirecionais em uma porta vertical considerando a
    interface z_i que separa a camada inferior (ar) da superior (fumaça).

    Equações baseadas no SFPE Handbook (3ª Ed., Eqs. 36.68–36.72).

    Parameters
    ----------
    largura : float
        Largura da porta [m].
    altura_porta : float
        Altura total da porta [m].
    altura_interface : float
        Altura da interface (Z_int) medida do piso [m].  0 ≤ z_i ≤ H.
    densidade_ar : float
        Densidade da camada inferior (ar) [kg/m³].
    densidade_fumaca : float
        Densidade da camada superior (fumaça) [kg/m³].
    C_d : float, opcional
        Coeficiente de descarga (padrão 0.68).
    m_dot_combustivel : float, opcional
        Vazão mássica de combustível (pirólise) [kg/s].

    Returns
    -------
    dict
        'z_n': altura do plano neutro [m]
        'm_dot_in': vazão mássica de ar frio entrando [kg/s]
        'm_dot_out': vazão mássica de fumaça saindo [kg/s]
    """
    # Se o solver testar uma densidade maior ou igual à do ar, não há empuxo térmico
    if densidade_fumaca >= densidade_ar:
         return {"z_n": altura_porta / 2.0, "m_dot_in": 0.0, "m_dot_out": 0.0}
    
    # Trava matemática suave para manter a interface dentro dos limites da porta
    if not (0 <= altura_interface <= altura_porta):
        altura_interface = max(0.0, min(altura_porta, altura_interface))

    g = 9.81
    W = largura
    H = altura_porta
    z_i = altura_interface
    drho = densidade_ar - densidade_fumaca
    rho_a = densidade_ar
    rho_f = densidade_fumaca

    # Coeficiente comum (⅔ C_d W √(2g))
    C0 = (2.0 / 3.0) * C_d * W * np.sqrt(2.0 * g)

    # ------------------------------------------------------------
    # Funções que calculam m_in e m_out para um dado z_n
    # ------------------------------------------------------------
    def fluxos(z_n):
        """Retorna (m_in, m_out) para um dado z_n."""
        if z_n < 0 or z_n > H:
            return np.inf, np.inf

        if z_n <= z_i:
            # CASO A: plano neutro abaixo da interface (ou na interface)
            # Entrada: ocorre de 0 a z_n na camada inferior (densidade = rho_a)
            # ΔP(z) = g [ drho*(H - z_i) - rho_a*(z_n - z) ]   para 0 ≤ z ≤ z_n
            # A integral ∫ √(ΔP) dz de 0 a z_n analítica:
            if z_n == 0:
                m_in = 0.0
            else:
                A = g * drho * (H - z_i)
                B = g * rho_a
                # termo = (2/(3*B)) * [ A^(3/2) - (A - B*z_n)^(3/2) ]
                # Cuidado: (A - B*z_n) pode ser negativo se z_n > A/B → nesse caso z_n estaria acima do ponto onde ΔP=0,
                # mas como z_n <= z_i, garante-se? Pode-se proteger com max(0, ...)
                A_val = max(0.0, A)
                inner = max(0.0, A_val - B * z_n)
                m_in = C0 * np.sqrt(rho_a) * (2.0 / (3.0 * B)) * (A_val ** 1.5 - inner ** 1.5)

            # Saída: ocorre de z_i a H na camada superior (densidade = rho_f)
            # ΔP(z) = g * drho * (H - z)   para z_i ≤ z ≤ H
            if z_i < H:
                m_out = C0 * np.sqrt(rho_f) * np.sqrt(g * drho) * (2.0 / 3.0) * (H - z_i) ** 1.5
            else:
                m_out = 0.0

        else:
            # CASO B: plano neutro acima da interface (z_n > z_i)
            # Entrada: ocorre de 0 a z_i (apenas na camada inferior)
            # ΔP(z) = g [ drho*(H - z_n) + rho_a*(z_n - z) ]   para 0 ≤ z ≤ z_i
            # (Eq. 36.71 do SFPE)
            if z_i == 0:
                m_in = 0.0
            else:
                A = g * drho * (H - z_n)
                B = g * rho_a
                # integral: (2/(3*B)) * [ (A + B*z_n)^(1.5) - (A + B*(z_n - z_i))^(1.5) ]
                termo1 = max(0.0, A + B * z_n)
                termo2 = max(0.0, A + B * (z_n - z_i))
                m_in = C0 * np.sqrt(rho_a) * (2.0 / (3.0 * B)) * (termo1 ** 1.5 - termo2 ** 1.5)

            # Saída: ocorre de z_n a H (na camada superior)
            # ΔP(z) = g * drho * (H - z)   para z_n ≤ z ≤ H
            if z_n < H:
                m_out = C0 * np.sqrt(rho_f) * np.sqrt(g * drho) * (2.0 / 3.0) * (H - z_n) ** 1.5
            else:
                m_out = 0.0

        return m_in, m_out

    # ------------------------------------------------------------
    # Resolver z_n pela condição m_in + m_dot_combustivel = m_out
    # ------------------------------------------------------------
    # Função objetivo f(z) = m_out - m_in - m_comb
    def f(z):
        m_in, m_out = fluxos(z)
        if np.isinf(m_in) or np.isinf(m_out):
            return np.inf
        return m_out - m_in - m_dot_combustivel

    # Chute inicial: usar solução sem interface (ou z_i como referência)
    # Aproveitamos a lógica da função original para dar um bom chute.
    if densidade_fumaca <= 0.0 or densidade_ar <= 0.0:
        return {"z_n": altura_porta / 2.0, "m_dot_in": 0.0, "m_dot_out": 0.0}
    r = densidade_fumaca / densidade_ar
    z0 = H * r ** (1 / 3) / (1 + r ** (1 / 3))
    # Se a interface estiver muito baixa, o plano neutro tende a cair; ajustamos o chute.
    if z_i < H:
        z0 = min(z0, z_i)  # nunca pode ser maior que a interface? Pode sim, mas é um chute.
    try:
        z_n_sol = newton(f, z0, tol=1e-8, maxiter=50)
        # Garante que fique no domínio
        z_n_sol = max(0.0, min(H, z_n_sol))
    except RuntimeError:
        # Fallback: usar solução analítica da função original (ignora interface)
        if densidade_ar <= 0.0 or densidade_fumaca <= 0.0:
            z_n_sol = altura_porta / 2.0
        else:
            z_n_sol = H / (1 + (densidade_ar / densidade_fumaca) ** (1.0 / 3.0))
    # Calcular fluxos com o z_n encontrado
    m_in_final, m_out_final = fluxos(z_n_sol)

    return {"z_n": z_n_sol, "m_dot_in": m_in_final, "m_dot_out": m_out_final}