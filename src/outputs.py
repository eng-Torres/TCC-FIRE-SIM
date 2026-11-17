"""
Módulo de Saída de Resultados
Gera gráficos no ecrã ou exporta os dados para ficheiro CSV.
Compatível com o modelo de parede de N nós (índice 5 = superfície).
"""
import numpy as np
import matplotlib.pyplot as plt
from src.fluidodinamica import calcular_fluxos_compartimento
from src.combustao import calcular_hrr_total_e_propriedades
from src.banco_construcao import MATERIAIS_CONSTRUCAO
from src.termica_parede import calcular_perda_calor_superior
from src.constantes import RHO_AMB, Y_O2_AMB, Y_CO2_AMB, Y_CO_AMB, T_AMB, SIGMA


def _calcular_campos(t, y, sala, C_P, R_AR, P_ATM):
    """
    Calcula as grandezas derivadas (temperatura, altura da interface,
    frações mássicas) a partir do estado y.
    y: [m_u, H_u, m_o2, m_co2, m_co, T_w0, ..., T_w(N-1)]
    """
    m_u, H_u, m_o2, m_co2, m_co = y[:5]
    T_w = y[5]

    T_u = H_u / np.maximum(m_u * C_P, 1e-12)   # <-- C_P incluído
    T_u_celsius = T_u - 273.15
    T_w_celsius = T_w - 273.15
    V_u = (m_u * R_AR * T_u) / P_ATM
    Z_int = sala["altura_z"] - (V_u / sala["area_piso"])
    Z_int = np.maximum(Z_int, 0.0)
    Y_o2  = m_o2  / np.maximum(m_u, 1e-12)
    Y_co2 = m_co2 / np.maximum(m_u, 1e-12)
    Y_co  = m_co  / np.maximum(m_u, 1e-12)
    return {
        't': t,
        'Z_int': Z_int,
        'T_u_celsius': T_u_celsius,
        'T_w_celsius': T_w_celsius,
        'Y_o2': Y_o2,
        'Y_co2': Y_co2,
        'Y_co': Y_co,
    }


def _calcular_fluxos(t, y, sala, aberturas, objetos_ativos,
                     C_P, R_AR, P_ATM, C_d=0.68):
    """Calcula m_dot_in e m_dot_out ao longo do tempo a partir do estado y."""
    n = len(t)
    m_dot_in_arr = np.zeros(n)
    m_dot_out_arr = np.zeros(n)

    rho_amb_local = P_ATM / (R_AR * T_AMB)
    nome_mat_parede = sala.get("material_parede", "gesso_16mm")
    props_parede = MATERIAIS_CONSTRUCAO.get(
        nome_mat_parede, MATERIAIS_CONSTRUCAO["gesso_16mm"]
    )
    geo_aberturas = [{
        "largura": vent["dimensoes"]["x"],
        "altura": vent["dimensoes"]["z"],
        "z_base": vent.get("posicao", {}).get("z", 0.0),
    } for vent in aberturas]

    for i in range(n):
        m_u = y[0, i] if y.ndim == 2 else y[0]
        H_u = y[1, i] if y.ndim == 2 else y[1]
        m_o2 = y[2, i] if y.ndim == 2 else y[2]
        T_w0 = y[5, i] if y.ndim == 2 else y[5]
        m_u = max(m_u, 0.001)
        m_o2 = max(m_o2, 0.0)
        T_u = max(200.0, min(H_u / (m_u * C_P), 2000.0))
        V_u = (m_u * R_AR * T_u) / P_ATM
        Z_int = sala["altura_z"] - V_u / sala["area_piso"]
        Z_int = max(Z_int, 0.0)
        rho_u = P_ATM / (R_AR * T_u)

        # 7 valores agora (incluindo D_fogo)
        Q_dot, m_dot_fuel, _, _, _, _, _ = \
            calcular_hrr_total_e_propriedades(t[i], objetos_ativos)

        # Limitação por O2 idêntica à integrada pelo solver (C_LOL)
        Y_O2_LOL = 0.01
        Y_O2_atual = m_o2 / m_u if m_u > 0 else 0.0
        C_LOL = 0.5 * (np.tanh(800.0 * (Y_O2_atual - Y_O2_LOL) - 4.0) + 1.0)
        m_dot_fuel *= C_LOL
        Q_c = 0.7 * Q_dot * C_LOL

        Q_loss, _ = calcular_perda_calor_superior(
            sala,
            Z_int,
            T_u,
            T_w0,
            props_parede["h_conv"],
            props_parede["emissividade"],
            sigma=SIGMA,
        )
        m_dot_exp_termica = (Q_c - Q_loss) / (C_P * T_AMB)

        fluxo = calcular_fluxos_compartimento(
            geo_aberturas,
            altura_interface=Z_int,
            densidade_ar=rho_amb_local,
            densidade_fumaca=rho_u,
            C_d=C_d,
            m_dot_combustivel=m_dot_fuel,
            m_dot_expansao=m_dot_exp_termica,
            razao_T_out=T_u / T_AMB,
        )

        m_dot_in_arr[i] = fluxo["m_dot_in"]
        m_dot_out_arr[i] = fluxo["m_dot_out"]

    return m_dot_in_arr, m_dot_out_arr

def plotar_graficos(t, y, sala, aberturas, C_P, R_AR, P_ATM, titulo="Simulação de Incêndio"):
    """Gera os quatro gráficos padrão no ecrã."""
    dados = _calcular_campos(t, y, sala, C_P, R_AR, P_ATM)
    Z_int = dados['Z_int']
    T_u_celsius = dados['T_u_celsius']
    T_w_celsius = dados['T_w_celsius']
    Y_o2 = dados['Y_o2']
    Y_co2 = dados['Y_co2']
    Y_co = dados['Y_co']

    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle(titulo, fontweight='bold')

    # Gráfico 1: Dinâmica da fumaça
    ax1.plot(t, Z_int, 'k', lw=2)
    ax1.axhline(sala["altura_z"], color='gray', ls='--', label='Teto')
    ax1.axhline(aberturas[0]["dimensoes"]["z"], color='red', ls=':', label='Altura da abertura')
    ax1.set_ylabel('Altura da interface (m)')
    ax1.set_title('Dinâmica da fumaça')
    ax1.legend(['Z_int', 'Teto', 'Altura da abertura'])

    # Gráfico 2: Temperaturas
    ax2.plot(t, T_u_celsius, 'firebrick', lw=2, label='Camada superior')
    ax2.plot(t, T_w_celsius, 'brown', lw=1.5, ls='--', label='Parede (superfície)')
    ax2.axhline(T_AMB-273.15, color='gray', ls=':', label=f'Ambiente ({T_AMB-273.15:.0f}°C)')
    ax2.axhline(600, color='darkorange', ls='--', label='Flashover 600°C')
    ax2.axhline(300, color='blue', ls=':', label='Ignição madeira 300°C')
    ax2.set_ylabel('Temperatura (°C)')
    ax2.set_title('Temperaturas do compartimento')
    ax2.legend()

    # Gráfico 3: Composição O₂ e CO₂
    ax3.plot(t, Y_o2*100, 'b', label='O₂')
    ax3.plot(t, Y_co2*100, 'g', label='CO₂')
    ax3.set_xlabel('Tempo (s)')
    ax3.set_ylabel('Fração mássica (%)')
    ax3.set_title('Composição na camada superior')
    ax3.legend()

    # Gráfico 4: CO
    ax4.plot(t, Y_co*1e6, 'purple', lw=2)  # ppm mássico
    ax4.set_xlabel('Tempo (s)')
    ax4.set_ylabel('CO (ppm)')
    ax4.set_title('Monóxido de Carbono')
    ax4.grid(True)

    plt.tight_layout()
    plt.show()

def exportar_dados(t, y, objetos_ativos, sala, aberturas, C_P, R_AR, P_ATM,
                   nome_ficheiro="resultados.csv", dt_print=10.0,
                   titulo="Simulação de Incêndio"):
    t_alvo = np.arange(0.0, t[-1] + dt_print, dt_print)

    # Interpola estado
    y_interp = np.array([np.interp(t_alvo, t, y[i]) for i in range(y.shape[0])])
    dados = _calcular_campos(t_alvo, y_interp, sala, C_P, R_AR, P_ATM)

    # --- Calcular vazões nos instantes de saída ---
    m_dot_in_arr, m_dot_out_arr = _calcular_fluxos(
    t_alvo, y_interp, sala, aberturas, objetos_ativos,
    C_P, R_AR, P_ATM, C_d=0.68)

    # Cabeçalho com as novas colunas
    cabecalho = (f"# Cenário: {titulo}\n"
                 "Tempo(s), Z_int(m), T_u(C), T_w(C), Y_O2, Y_CO2, Y_CO, m_dot_in(kg/s), m_dot_out(kg/s)")
    linhas = []
    for i in range(len(t_alvo)):
        linhas.append(
            f"{t_alvo[i]:.2f}, {dados['Z_int'][i]:.4f}, {dados['T_u_celsius'][i]:.2f}, "
            f"{dados['T_w_celsius'][i]:.2f}, {dados['Y_o2'][i]:.6f}, "
            f"{dados['Y_co2'][i]:.6f}, {dados['Y_co'][i]:.6f}, "
            f"{m_dot_in_arr[i]:.6f}, {m_dot_out_arr[i]:.6f}"
        )

    with open(nome_ficheiro, 'w') as f:
        f.write(cabecalho + "\n")
        f.write("\n".join(linhas))
    print(f"Dados exportados para {nome_ficheiro} (dt = {dt_print} s)")