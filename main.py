"""
Script principal para simulação e visualização do modelo de duas zonas.
O cenário é selecionado pela variável CENARIO.
"""
from src.geometria import carregar_geometria_padrao, carregar_geometria_teste2_1991
from src.edos_solver import simular_incendio_sequencial
from src.constantes import C_P, R_AR, P_ATM
from src.outputs import plotar_graficos, exportar_dados

# ============================================================
# Escolha do cenário (pode vir de input() se preferir)
# ============================================================
CENARIO = "teste2_1991"   # opções: "teste2_1991", "quarto_padrao"
MODO_SAIDA = 'plot'       # 'plot' ou 'csv'

def configurar_cenario(nome: str):
    """
    Retorna um dicionário com todos os parâmetros do cenário escolhido.
    Chaves:
        sala, aberturas, objeto_inicial, outros_materiais,
        objetos_geometria, titulo, t_final
    """
    if nome == "teste2_1991":
        sala, _, aberturas = carregar_geometria_teste2_1991()
        objeto_inicial = {"nome_material": "teste2_hrr_medida"}
        outros_materiais = []
        objetos_geometria = [{"nome_material": "teste2_hrr_medida", "area": 1.5}]
        titulo = "Teste 2 – Peacock et al. (1991)"
        t_final = 2070
    elif nome == "quarto_padrao":
        sala, objetos, aberturas = carregar_geometria_padrao()
        cama = next(obj for obj in objetos if obj["nome"] == "Cama de Casal")
        guarda_roupa = next(obj for obj in objetos if obj["nome"] == "Guarda-Roupa")
        objeto_inicial = {"nome_material": cama["tipo_material"]}
        outros_materiais = [
            {"nome_material": guarda_roupa["tipo_material"], "T_ign": guarda_roupa["T_ign"]}
        ]
        objetos_geometria = [
            {"nome_material": cama["tipo_material"], "area": 2.8},
            {"nome_material": guarda_roupa["tipo_material"], "area": 0.9}
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
    modelo_pluma = 'zukoski'  # 'zukoski', 'mccaffrey', 'heskestad'

    cfg = configurar_cenario(CENARIO)
    print(f"Cenário: {CENARIO}  |  Pluma: {modelo_pluma}")
    print(f"Iniciando simulação de {cfg['t_final']} segundos. Aguarde...")

    t, y = simular_incendio_sequencial(
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
        exportar_dados(t, y, cfg['sala'], cfg['aberturas'],
                       C_P, R_AR, P_ATM,
                       nome_ficheiro=nome_ficheiro, dt_print=10.0,
                       titulo=cfg['titulo'])
    else:
        print("Modo de saída desconhecido. Use 'plot' ou 'csv'.")


if __name__ == "__main__":
    
    rodar_teste()