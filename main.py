"""
Script principal do simulador de incêndio (modelo de duas zonas).

Uso:
    python main.py --listar
    python main.py --listar-materiais
    python main.py --cenario teste1_Babrauskas
    python main.py --cenario quarto_padrao --pluma heskestad
    python main.py --cenario teste1_Babrauskas --t-final 1200 --saida ambos

Cenários ficam em ./cenarios/*.yaml (ou .json) e podem ser versionados no Git.
"""
import argparse
import sys

from src.io_cenario import (
    carregar_cenario, listar_cenarios, listar_materiais,
    localizar_cenario, registrar_materiais_customizados,
    CENARIOS_DIR, ErroCenario,
)
from src.edos_solver import simular_incendio_sequencial
from src.constantes import C_P, R_AR, P_ATM
from src.outputs import plotar_graficos, exportar_dados


def _montar_objetos_para_solver(objetos_json):
    """
    Converte a lista unificada de objetos do cenário para a interface
    esperada por simular_incendio_sequencial.
    """
    objeto_inicial = None
    outros = []
    geometria = []

    for obj in objetos_json:
        geometria.append({
            "nome_material": obj["nome_material"],
            "area": obj["area"],
        })
        if obj["is_initial"]:
            objeto_inicial = {"nome_material": obj["nome_material"]}
        else:
            outros.append({
                "nome_material": obj["nome_material"],
                "T_ign": obj["T_ign"],
            })

    return objeto_inicial, outros, geometria


def rodar(cfg, args):
    sim = cfg["simulacao"]

    # Overrides via CLI
    if args.pluma:
        sim["modelo_pluma"] = args.pluma
    if args.t_final is not None:
        sim["t_final"] = args.t_final
    if args.saida:
        cfg["saida"]["modo"] = args.saida
    if args.nome_saida:
        cfg["saida"]["nome_ficheiro"] = args.nome_saida

    print(f"Cenário: {cfg['nome']}")
    if cfg["descricao"]:
        print(f"  Descrição: {cfg['descricao'].strip()}")
    if cfg["referencia"]:
        print(f"  Referência: {cfg['referencia'].strip()}")
    print(f"  Pluma: {sim['modelo_pluma']}  |  t_final: {sim['t_final']} s")
    print(f"  N_NOS={sim['n_nos_parede']}  MAX_STEP={sim['max_step']}  "
          f"PPS={sim['pontos_por_segundo']}")
    print(f"  Material da parede: {cfg['sala']['material_parede']}")
    print("Iniciando simulação...")

    objeto_inicial, outros, geometria = _montar_objetos_para_solver(cfg["objetos"])

    t, y, objetos_ativos = simular_incendio_sequencial(
        cfg["sala"],
        cfg["aberturas"],
        objeto_inicial=objeto_inicial,
        outros_objetos=outros,
        objetos_geometria=geometria,
        t_final=sim["t_final"],
        modelo_pluma=sim["modelo_pluma"],
        max_step=sim["max_step"],
        pontos_por_segundo=sim["pontos_por_segundo"],
        n_nos=sim["n_nos_parede"],
    )

    print("Simulação concluída! Processando resultados...")

    modo = cfg["saida"]["modo"]
    titulo = f"{cfg['nome']} — pluma {sim['modelo_pluma']}"

    if modo in ("plot", "ambos"):
        plotar_graficos(t, y, cfg["sala"], cfg["aberturas"],
                        C_P, R_AR, P_ATM, titulo=titulo)

    if modo in ("csv", "ambos"):
        exportar_dados(t, y, objetos_ativos, cfg["sala"], cfg["aberturas"],
                       C_P, R_AR, P_ATM,
                       nome_ficheiro=cfg["saida"]["nome_ficheiro"],
                       dt_print=sim["dt_print"],
                       titulo=titulo)


def main():
    parser = argparse.ArgumentParser(
        description="Simulador de incêndio — modelo de duas zonas"
    )
    parser.add_argument(
        "--cenario", "-c",
        help="Nome do cenário (arquivo em ./cenarios/ sem extensão)"
    )
    parser.add_argument(
        "--listar", action="store_true",
        help="Lista cenários disponíveis e sai"
    )
    parser.add_argument(
        "--listar-materiais", action="store_true", dest="listar_materiais",
        help="Lista materiais embutidos (banco) e sai"
    )
    parser.add_argument(
        "--pluma", choices=["heskestad", "zukoski", "mccaffrey"],
        help="Override do modelo de pluma"
    )
    parser.add_argument(
        "--t-final", type=float, dest="t_final",
        help="Override do tempo final [s]"
    )
    parser.add_argument(
        "--saida", choices=["csv", "plot", "ambos"],
        help="Override do modo de saída"
    )
    parser.add_argument(
        "--nome-saida", dest="nome_saida",
        help="Override do nome do arquivo CSV"
    )

    args = parser.parse_args()

    if args.listar:
        cenarios = listar_cenarios()
        if not cenarios:
            print(f"Nenhum cenário encontrado em {CENARIOS_DIR}")
        else:
            print("Cenários disponíveis:")
            for c in cenarios:
                print(f"  - {c}")
        sys.exit(0)

    if args.listar_materiais:
        m = listar_materiais()
        print("Materiais combustíveis (banco embutido):")
        for n in m["combustiveis"]:
            print(f"  - {n}")
        print("Materiais de construção (banco embutido):")
        for n in m["construcao"]:
            print(f"  - {n}")
        sys.exit(0)

    if not args.cenario:
        parser.error("--cenario é obrigatório (ou use --listar)")

    try:
        caminho = localizar_cenario(args.cenario)
        cfg = carregar_cenario(caminho)
    except ErroCenario as e:
        print(f"[ERRO] {e}", file=sys.stderr)
        sys.exit(1)

    registrar_materiais_customizados(cfg)
    rodar(cfg, args)


if __name__ == "__main__":
    main()
