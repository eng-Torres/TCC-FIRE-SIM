"""
Módulo de Saída de Resultados
Gera gráficos no ecrã ou exporta os dados para ficheiro CSV.
Compatível com o modelo de parede de N nós (índice 5 = superfície).
"""
import numpy as np
import matplotlib.pyplot as plt

def _calcular_campos(t, y, sala, C_P, R_AR, P_ATM):
    """
    Calcula as grandezas derivadas (temperatura, altura da interface,
    frações mássicas) a partir do estado y.
    y: [m_u, H_u, m_o2, m_co2, m_co, T_w0, ..., T_w(N-1)]
    """
    m_u, H_u, m_o2, m_co2, m_co = y[:5]
    T_w = y[5]                       # temperatura da superfície (primeiro nó)

    T_u = H_u / np.maximum(m_u, 1e-6)
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
    ax2.axhline(20, color='gray', ls=':', label='Ambiente (20°C)')
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

def exportar_dados(t, y, sala, aberturas, C_P, R_AR, P_ATM,
                   nome_ficheiro="resultados.csv", dt_print=10.0,
                   titulo="Simulação de Incêndio"):
    """
    Escreve as variáveis de estado num ficheiro CSV, amostradas de dt_print em dt_print segundos.
    O parâmetro 'titulo' é incluído como comentário no cabeçalho, mas não afeta os dados.
    """
    t_alvo = np.arange(0.0, t[-1] + dt_print, dt_print)

    # Interpola os valores para os instantes desejados
    y_interp = np.array([np.interp(t_alvo, t, y[i]) for i in range(y.shape[0])])
    dados = _calcular_campos(t_alvo, y_interp, sala, C_P, R_AR, P_ATM)

    cabecalho = f"# Cenário: {titulo}\nTempo(s), Z_int(m), T_u(C), T_w(C), Y_O2, Y_CO2, Y_CO"
    linhas = []
    for i in range(len(t_alvo)):
        linhas.append(
            f"{t_alvo[i]:.2f}, {dados['Z_int'][i]:.4f}, {dados['T_u_celsius'][i]:.2f}, "
            f"{dados['T_w_celsius'][i]:.2f}, {dados['Y_o2'][i]:.6f}, "
            f"{dados['Y_co2'][i]:.6f}, {dados['Y_co'][i]:.6f}"
        )

    with open(nome_ficheiro, 'w') as f:
        f.write(cabecalho + "\n")
        f.write("\n".join(linhas))
    print(f"Dados exportados para {nome_ficheiro} (dt = {dt_print} s)")