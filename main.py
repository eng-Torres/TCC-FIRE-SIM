import numpy as np
import matplotlib.pyplot as plt

# Importações baseadas na sua arquitetura modular
from src.geometria import carregar_geometria_padrao
from src.edos_solver import simular_incendio

def rodar_teste():
    print("Carregando cenário padrão...")
    sala, objetos, aberturas = carregar_geometria_padrao()
    
    tempo_simulacao = 600 # 10 minutos
    print(f"Iniciando simulação de {tempo_simulacao} segundos. Por favor, aguarde...")
    
    # Roda o motor matemático
    solucao = simular_incendio(sala, aberturas, crescimento="medio", t_final=tempo_simulacao)
    
    if not solucao.success:
        print("A simulação falhou:", solucao.message)
        return
        
    print("Simulação concluída com sucesso! Processando resultados...")
    
    # Extrair os arrays de tempo e variáveis de estado
    tempos = solucao.t
    m_u = solucao.y[0]
    H_u = solucao.y[1]
    
    # Constantes
    C_P = 1.005
    R_AR = 0.287
    P_ATM = 101.325
    
    # Recalcular variáveis dependentes para o gráfico
    T_u = H_u / (m_u * C_P)
    T_u_celsius = T_u - 273.15 # Conversão para Celsius para leitura intuitiva
    
    V_u = (m_u * R_AR * T_u) / P_ATM
    Z_int = sala["altura_z"] - (V_u / sala["area_piso"])
    Z_int = np.maximum(Z_int, 0.0) # Clamping visual para não passar do chão
    
    # --- PLOTAGEM DOS RESULTADOS ---
    plt.style.use('seaborn-v0_8-darkgrid')
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Gráfico 1: Altura da Camada de Fumaça
    ax1.plot(tempos, Z_int, color='black', linewidth=2, label='Interface da Fumaça ($Z_{int}$)')
    ax1.axhline(y=sala["altura_z"], color='gray', linestyle='--', label='Teto')
    ax1.axhline(y=aberturas[0]["dimensoes"]["z"], color='red', linestyle=':', label='Altura da Porta')
    ax1.set_ylabel('Altura (m)', fontsize=12)
    ax1.set_title('Dinâmica da Fumaça no Compartimento', fontsize=14, fontweight='bold')
    ax1.set_ylim(0, sala["altura_z"] + 0.2)
    ax1.legend(loc='lower left')
    
    # Gráfico 2: Temperatura da Camada de Fumaça
    ax2.plot(tempos, T_u_celsius, color='firebrick', linewidth=2, label='Temperatura da Fumaça ($T_u$)')
    ax2.axhline(y=600, color='darkorange', linestyle='--', label='Risco de Flashover (600°C)')
    ax2.set_xlabel('Tempo (s)', fontsize=12)
    ax2.set_ylabel('Temperatura (°C)', fontsize=12)
    ax2.legend(loc='upper left')
    
    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    rodar_teste()