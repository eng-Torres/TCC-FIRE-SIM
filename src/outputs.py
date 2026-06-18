"""
Módulo de Saída de Resultados
Gera gráficos no ecrã ou exporta os dados para ficheiro CSV.
"""
import numpy as np
import matplotlib.pyplot as plt

def _extrair_campos(t, y, sala, aberturas, C_P, R_AR, P_ATM):
    """Extrai as variáveis comuns a todas as saídas."""
    m_u, H_u, m_o2, m_co2, m_co, T_w, D_fogo = y
    T_u = H_u / (m_u * C_P)
    T_u_celsius = T_u - 273.15
    T_w_celsius = T_w - 273.15
    V_u = (m_u * R_AR * T_u) / P_ATM
    Z_int = sala["altura_z"] - (V_u / sala["area_piso"])
    Z_int = np.maximum(Z_int, 0.0)
    Y_o2  = m_o2  / m_u
    Y_co2 = m_co2 / m_u
    Y_co  = m_co  / m_u
    return {
        't': t,
        'Z_int': Z_int,
        'T_u_celsius': T_u_celsius,
        'T_w_celsius': T_w_celsius,
        'Y_o2': Y_o2,
        'Y_co2': Y_co2,
        'Y_co': Y_co,
        'D_fogo': D_fogo,
        'sala': sala,
        'aberturas': aberturas
    }

def plotar_graficos(t, y, sala, aberturas, C_P, R_AR, P_ATM):
    """
    Gera os quatro gráficos padrão no ecrã.
    """
    dados = _extrair_campos(t, y, sala, aberturas, C_P, R_AR, P_ATM)
    Z_int = dados['Z_int']
    T_u_celsius = dados['T_u_celsius']
    T_w_celsius = dados['T_w_celsius']
    Y_o2 = dados['Y_o2']
    Y_co2 = dados['Y_co2']
    Y_co = dados['Y_co']
    sala = dados['sala']
    aberturas = dados['aberturas']

    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('Incêndio com ignição sequencial (Cama → Guarda-Roupa)', fontweight='bold')

    # Gráfico 1: Dinâmica da fumaça
    ax1.plot(t, Z_int, 'k', lw=2)
    ax1.axhline(sala["altura_z"], color='gray', ls='--', label='Teto')
    ax1.axhline(aberturas[0]["dimensoes"]["z"], color='red', ls=':', label='Altura da porta')
    ax1.set_ylabel('Altura da interface (m)')
    ax1.set_title('Dinâmica da fumaça')
    ax1.legend(['Z_int', 'Teto', 'Altura da porta'])

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
                   nome_ficheiro="resultados.csv", dt_print=10.0):
    """
    Escreve as variáveis de estado num ficheiro CSV, amostradas de dt_print em dt_print segundos.
    """
    # Gerar os instantes de saída (desde t=0 até ao final, com passo dt_print)
    t_alvo = np.arange(0.0, t[-1] + dt_print, dt_print)
    
    # Interpolar cada variável para os instantes alvo
    m_u    = np.interp(t_alvo, t, y[0])
    H_u    = np.interp(t_alvo, t, y[1])
    m_o2   = np.interp(t_alvo, t, y[2])
    m_co2  = np.interp(t_alvo, t, y[3])
    m_co   = np.interp(t_alvo, t, y[4])
    T_w    = np.interp(t_alvo, t, y[5])
    D_fogo = np.interp(t_alvo, t, y[6])
    
    # Calcular as grandezas derivadas
    T_u = H_u / (m_u * C_P)
    T_u_celsius = T_u - 273.15
    T_w_celsius = T_w - 273.15
    V_u = (m_u * R_AR * T_u) / P_ATM
    Z_int = sala["altura_z"] - (V_u / sala["area_piso"])
    Z_int = np.maximum(Z_int, 0.0)
    Y_o2  = m_o2  / m_u
    Y_co2 = m_co2 / m_u
    Y_co  = m_co  / m_u
    
    # Escrever o CSV
    cabecalho = ("Tempo(s), Z_int(m), T_u(C), T_w(C), Y_O2, Y_CO2, Y_CO, D_fogo(m)")
    linhas = []
    for i in range(len(t_alvo)):
        linhas.append(f"{t_alvo[i]:.2f}, {Z_int[i]:.4f}, {T_u_celsius[i]:.2f}, "
                      f"{T_w_celsius[i]:.2f}, {Y_o2[i]:.6f}, "
                      f"{Y_co2[i]:.6f}, {Y_co[i]:.6f}, {D_fogo[i]:.4f}")
    
    with open(nome_ficheiro, 'w') as f:
        f.write(cabecalho + "\n")
        f.write("\n".join(linhas))
    print(f"Dados exportados para {nome_ficheiro} (dt = {dt_print} s)")