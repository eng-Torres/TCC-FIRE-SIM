"""
Script principal para simulação e visualização do modelo de duas zonas.
"""
from src.geometria import carregar_geometria_padrao, carregar_geometria_teste1_Babrauskas
from src.edos_solver import simular_incendio_sequencial
from src.constantes import C_P, R_AR, P_ATM
from src.outputs import plotar_graficos, exportar_dados

CENARIO = "teste1_Babrauskas"   # opções: "teste1_Babrauskas", "quarto_padrao"
MODO_SAIDA = 'csv' # plot , csv


def configurar_cenario(nome: str):
    if nome == "teste1_Babrauskas":
        sala, _, aberturas = carregar_geometria_teste1_Babrauskas()
        objeto_inicial = {"nome_material": "teste1_hrr_medida"}
        outros_materiais = []
        # Área real do loveseat F31 = 1.5 × 0.8 = 1.2 m²
        objetos_geometria = [{"nome_material": "teste1_hrr_medida", "area": 1.2}]
        titulo = "Teste 1 – Babrauskas (1984)"
        t_final = 2070

    elif nome == "quarto_padrao":
        sala, objetos, aberturas = carregar_geometria_padrao()
        cama = next(obj for obj in objetos if obj["nome"] == "Cama de Casal")
        guarda_roupa = next(obj for obj in objetos if obj["nome"] == "Guarda-Roupa")
        tapete_de_papelao = next(obj for obj in objetos if obj["nome"] == "Tapete de Papelao")
        cadeira = next(obj for obj in objetos if obj["nome"] == "cadeira")
        
        objeto_inicial = {"nome_material": cama["tipo_material"]}
        outros_materiais = [
            {"nome_material": guarda_roupa["tipo_material"], "T_ign": guarda_roupa["T_ign"]},
            {"nome_material": tapete_de_papelao["tipo_material"], "T_ign": tapete_de_papelao["T_ign"]},
            {"nome_material": cadeira["tipo_material"], "T_ign": cadeira["T_ign"]}
        ]
        objetos_geometria = [
            {"nome_material": cama["tipo_material"], "area": 2.8},
            {"nome_material": guarda_roupa["tipo_material"], "area": 0.9},
            {"nome_material": tapete_de_papelao["tipo_material"], "area": 2.0},
            {"nome_material": cadeira["tipo_material"], "area": 0.25}
        ]
        titulo = "Incêndio no Quarto Padrão – Ignição Sequencial"
        t_final = 600
    else:
        raise ValueError(f"Cenário '{nome}' não reconhecido.")

    return {
        'sala': sala,
        'aberturas': aberturas,
        'objeto_inicial': objeto_inicial,
        'outros_materiais': outros_materiais,
        'objetos_geometria': objetos_geometria,
        'titulo': titulo,
        't_final': t_final
    }


def rodar_teste():
    modelo_pluma = 'mccaffrey' # zukoski , heskestad , mccaffrey

    cfg = configurar_cenario(CENARIO)
    print(f"Cenário: {CENARIO}  |  Pluma: {modelo_pluma}")
    print(f"Iniciando simulação de {cfg['t_final']} segundos. Aguarde...")

    t, y, objetos_ativos = simular_incendio_sequencial(
        cfg['sala'], cfg['aberturas'],
        objeto_inicial=cfg['objeto_inicial'],
        outros_objetos=cfg['outros_materiais'],
        objetos_geometria=cfg['objetos_geometria'],
        t_final=cfg['t_final'],
        modelo_pluma=modelo_pluma,
        crescimento_t2="medio"
    )

    print("Simulação concluída! Processando resultados...")

    if MODO_SAIDA == 'plot':
        plotar_graficos(t, y, cfg['sala'], cfg['aberturas'],
                        C_P, R_AR, P_ATM, titulo=cfg['titulo'])
    elif MODO_SAIDA == 'csv':
        nome_ficheiro = f"resultados_{CENARIO}_{modelo_pluma}.csv"
        exportar_dados(t, y, objetos_ativos, cfg['sala'], cfg['aberturas'],
                       C_P, R_AR, P_ATM,
                       nome_ficheiro=nome_ficheiro, dt_print=10.0,
                       titulo=cfg['titulo'])
    else:
        print("Modo de saída desconhecido. Use 'plot' ou 'csv'.")


if __name__ == "__main__":
    rodar_teste()