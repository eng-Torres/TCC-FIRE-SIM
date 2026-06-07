"""
Script principal para simulação e visualização do modelo de duas zonas.
Agora com ignição sequencial de objetos (cama → guarda-roupa).
"""

import numpy as np
import matplotlib.pyplot as plt
from src.geometria import carregar_geometria_padrao
from src.edos_solver import simular_incendio_sequencial

def rodar_teste():
    print("Carregando cenário padrão...")
    sala, objetos, aberturas = carregar_geometria_padrao()

    # Identificar o objeto inicial (Cama de Casal) e os demais
    # A cama é o objeto com tipo_material = 'espuma_poliuretano'
    cama = next(obj for obj in objetos if obj["nome"] == "Cama de Casal")
    guarda_roupa = next(obj for obj in objetos if obj["nome"] == "Guarda-Roupa")

    objeto_inicial = {"nome_material": cama["tipo_material"]}
    outros_materiais = [guarda_roupa["tipo_material"]]

    print(f"Foco inicial: {cama['nome']} ({cama['tipo_material']})")
    print(f"Objeto secundário: {guarda_roupa['nome']} ({guarda_roupa['tipo_material']})")

    tempo_simulacao = 600  # 10 minutos
    print(f"Iniciando simulação de {tempo_simulacao} segundos. Aguarde...")

    t, y = simular_incendio_sequencial(
        sala, aberturas,
        objeto_inicial=objeto_inicial,
        outros_objetos=outros_materiais,
        t_final=tempo_simulacao,
        modelo_pluma='zukoski', # padrão é zukoski, e a outra opção é heskestad (NOT FULLY FUNCIONAL, NEED FIX)
        crescimento_t2="medio"  # fallback se material não existir
    )

    print("Simulação concluída! Processando resultados...")

    m_u   = y[0]
    H_u   = y[1]
    m_o2  = y[2]
    m_co2 = y[3]
    m_co  = y[4]
    T_w   = y[5]

    C_P = 1.005
    R_AR = 0.287
    P_ATM = 101.325

    T_u = H_u / (m_u * C_P)
    T_u_celsius = T_u - 273.15
    T_w_celsius = T_w -273.15
    V_u = (m_u * R_AR * T_u) / P_ATM
    Z_int = sala["altura_z"] - (V_u / sala["area_piso"])
    Z_int = np.maximum(Z_int, 0.0)

    Y_o2  = m_o2  / m_u
    Y_co2 = m_co2 / m_u
    Y_co  = m_co  / m_u
    

    # Gráficos
    plt.style.use('seaborn-v0_8-darkgrid')
    fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, figsize=(12, 8))
    fig.suptitle('Incêndio com ignição sequencial (Cama → Guarda-Roupa)', fontweight='bold')

    ax1.plot(t, Z_int, 'k', lw=2)
    ax1.axhline(sala["altura_z"], color='gray', ls='--')
    ax1.axhline(aberturas[0]["dimensoes"]["z"], color='red', ls=':')
    ax1.set_ylabel('Altura da interface (m)')
    ax1.set_title('Dinâmica da fumaça')
    ax1.legend(['Z_int', 'Teto', 'Altura da porta'])

    ax2.plot(t, T_u_celsius, 'firebrick', lw=2, label='Camada superior')
    ax2.plot(t, T_w_celsius, 'brown', lw=1.5, ls='--', label='Parede (superfície)')
    ax2.axhline(20, color='gray', ls=':', label='Ambiente (20°C)')
    ax2.axhline(600, color='darkorange', ls='--', label='Flashover 600°C')
    ax2.axhline(300, color='blue', ls=':', label='Ignição madeira 300°C')
    ax2.set_ylabel('Temperatura (°C)')
    ax2.set_title('Temperaturas do compartimento')
    ax2.legend()

    ax3.plot(t, Y_o2*100, 'b', label='O₂')
    ax3.plot(t, Y_co2*100, 'g', label='CO₂')
    ax3.set_xlabel('Tempo (s)')
    ax3.set_ylabel('Fração mássica (%)')
    ax3.set_title('Composição na camada superior')
    ax3.legend()

    ax4.plot(t, Y_co*1e6, 'purple', lw=2)  # ppm mássico
    ax4.set_xlabel('Tempo (s)')
    ax4.set_ylabel('CO (ppm)')
    ax4.set_title('Monóxido de Carbono')
    ax4.grid(True)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    rodar_teste()