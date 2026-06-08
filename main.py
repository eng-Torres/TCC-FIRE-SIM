"""
Script principal para simulação e visualização do modelo de duas zonas.
"""
import numpy as np
import matplotlib.pyplot as plt
from src.geometria import carregar_geometria_padrao
from src.edos_solver import simular_incendio_sequencial
from src.constantes import C_P, R_AR, P_ATM, T_AMB
from src.outputs import plotar_graficos, exportar_dados

# Escolhe o modo de saída: 'plot' para gráficos, 'csv' para ficheiro de texto
MODO_SAIDA = 'plot'   # alterna para 'csv' quando quiseres exportar os dados

def rodar_teste():
    print("Carregando cenário padrão...")
    sala, objetos, aberturas = carregar_geometria_padrao()

    # Identificar objetos
    cama = next(obj for obj in objetos if obj["nome"] == "Cama de Casal")
    guarda_roupa = next(obj for obj in objetos if obj["nome"] == "Guarda-Roupa")
    tapete = next(obj for obj in objetos if obj["nome"] == "Tapete de Papelao")

    objeto_inicial = {"nome_material": cama["tipo_material"]}
    outros_materiais = [
        {"nome_material": tapete["tipo_material"], "T_ign": tapete["T_ign"]},
        {"nome_material": guarda_roupa["tipo_material"], "T_ign": guarda_roupa["T_ign"]}
    ]

    # Construir geometria dos objetos (áreas)
    objetos_geometria = []
    for obj in objetos:
        area = obj["dimensoes"]["x"] * obj["dimensoes"]["y"]
        objetos_geometria.append({"nome_material": obj["tipo_material"], "area": area})

    print(f"Foco inicial: {cama['nome']} ({cama['tipo_material']})")
    tempo_simulacao = 600
    print(f"Iniciando simulação de {tempo_simulacao} segundos. Aguarde...")

    t, y = simular_incendio_sequencial(
        sala, aberturas,
        objeto_inicial=objeto_inicial,
        outros_objetos=outros_materiais,
        objetos_geometria=objetos_geometria,
        t_final=tempo_simulacao,
        modelo_pluma='heskestad',
        crescimento_t2="medio"
    )

    print("Simulação concluída! Processando resultados...")

    if MODO_SAIDA == 'plot':
        plotar_graficos(t, y, sala, aberturas, C_P, R_AR, P_ATM)
    elif MODO_SAIDA == 'csv':
        exportar_dados(t, y, sala, aberturas, C_P, R_AR, P_ATM, nome_ficheiro="resultados_simulacao.csv")
    else:
        print("Modo de saída desconhecido. Use 'plot' ou 'csv'.")

if __name__ == "__main__":
    rodar_teste()