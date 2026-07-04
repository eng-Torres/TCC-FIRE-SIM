"""
Script principal para simulação e visualização do modelo de duas zonas.
"""
import numpy as np
import matplotlib.pyplot as plt
from src.geometria import carregar_geometria_padrao
from src.geometria import carregar_geometria_teste2_1991
from src.edos_solver import simular_incendio_sequencial
from src.constantes import C_P, R_AR, P_ATM, T_AMB
from src.outputs import plotar_graficos, exportar_dados

# Escolhe o modo de saída: 'plot' para gráficos, 'csv' para ficheiro de texto
MODO_SAIDA = 'plot'   # alterna para 'csv' quando quiseres exportar os dados

def rodar_teste():

    modelo_pluma='zukoski'  # modelos de pluma mccaffrey,heskestad,zukoski

    print("Carregando cenário padrão...")
    sala, objetos, aberturas = carregar_geometria_teste2_1991()
    # Identificar objetos
    cesto = next(obj for obj in objetos if obj["nome"] == "Cesto de Lixo")
    loveseat = next(obj for obj in objetos if obj["nome"] == "Love Seat")

    objeto_inicial = {"nome_material": cesto["tipo_material"]}
    outros_materiais = [
        {"nome_material": loveseat["tipo_material"], "T_ign": loveseat["T_ign"]}
    ]

    # Construir geometria dos objetos (áreas)
    objeto_inicial = {"nome_material": "teste2_hrr_medida"}
    outros_materiais = []   # sem objetos secundários

    objetos_geometria = [
        {"nome_material": "teste2_hrr_medida", "area": 1.5}  
    ]
    tempo_simulacao = 2070
    print(f"Iniciando simulação de {tempo_simulacao} segundos. Aguarde...")

    t, y = simular_incendio_sequencial(
        sala, aberturas,
        objeto_inicial=objeto_inicial,
        outros_objetos=outros_materiais,
        objetos_geometria=objetos_geometria,
        t_final=tempo_simulacao,
        modelo_pluma=modelo_pluma, 
        crescimento_t2="medio"
    )

    print("Simulação concluída! Processando resultados...")

    if MODO_SAIDA == 'plot':
        plotar_graficos(t, y, sala, aberturas, C_P, R_AR, P_ATM)
    elif MODO_SAIDA == 'csv':
        nome_ficheiro = f"resultados_{modelo_pluma}.csv"
        exportar_dados(t, y, sala, aberturas, C_P, R_AR, P_ATM, nome_ficheiro=nome_ficheiro,dt_print=10.0)
    else:
        print("Modo de saída desconhecido. Use 'plot' ou 'csv'.")

if __name__ == "__main__":
    rodar_teste()