"""
Validação do simulador contra dados experimentais e/ou CFAST.

Uso:
    python valid.py --cenario teste1_Babrauskas
    python valid.py --cenario teste1_Babrauskas --t-max 1200

O cenário YAML precisa ter uma seção 'validacao' (ver exemplo).
Todos os caminhos em 'validacao' são relativos a ./cenarios/.
"""
import argparse
import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

from src.io_cenario import carregar_cenario, localizar_cenario, ErroCenario
from src.io_dados import ler_serie, interpolar_em_grid

CENARIOS_DIR = Path(__file__).resolve().parent / "cenarios"


# ---------------------------------------------------------------------------
# Métricas de erro
# ---------------------------------------------------------------------------
def metricas(y_sim, y_ref, nome):
    """RMSE, MAE, Max, viés — em unidades físicas. Sem MAPE (instável)."""
    mask = np.isfinite(y_sim) & np.isfinite(y_ref)
    if mask.sum() == 0:
        return None
    d = y_sim[mask] - y_ref[mask]
    return {
        'nome': nome,
        'n':    int(mask.sum()),
        'rmse': float(np.sqrt(np.mean(d ** 2))),
        'mae':  float(np.mean(np.abs(d))),
        'max':  float(np.max(np.abs(d))),
        'bias': float(np.mean(d)),
    }


def imprimir_metricas(tabela):
    print(f"\n{'Série':<35} {'n':>5} {'RMSE':>9} {'MAE':>9} {'Max':>9} {'Viés':>9}")
    print('-' * 80)
    for m in tabela:
        print(f"{m['nome']:<35} {m['n']:>5d} "
              f"{m['rmse']:>9.3f} {m['mae']:>9.3f} "
              f"{m['max']:>9.3f} {m['bias']:>9.3f}")


# ---------------------------------------------------------------------------
# Leitura dos arquivos listados no YAML
# ---------------------------------------------------------------------------
def _resolver_caminho(rel):
    p = Path(rel)
    if not p.is_absolute():
        p = CENARIOS_DIR / p
    return p


def carregar_series(cfg):
    """
    Lê experimental, CFAST e modelos a partir da seção 'validacao'.

    Arquivos ausentes são tratados como fontes opcionais: um aviso é
    impresso e a fonte é ignorada. O programa só falha se:
      - não houver seção 'validacao';
      - não sobrar nenhuma série válida (nem experimental, nem modelos).
    """
    val = cfg.get('validacao')
    if not val:
        raise ErroCenario(
            "Cenário não tem seção 'validacao'. Adicione-a para usar valid.py."
        )

    resultado = {'experimental': None, 'cfast': [], 'modelos': []}

    # --- Experimental (obrigatório) ---
    exp = val.get('experimental')
    if exp is None:
        raise ErroCenario("Seção 'validacao.experimental' é obrigatória.")

    caminho_exp = _resolver_caminho(exp['arquivo'])
    if not caminho_exp.is_file():
        raise ErroCenario(
            f"Arquivo experimental não encontrado: {caminho_exp}\n"
            f"Coloque o arquivo em ./cenarios/ ou ajuste o caminho no YAML."
        )

    try:
        serie = ler_serie(caminho_exp, formato=exp.get('formato', 'auto'))
    except Exception as e:
        raise ErroCenario(
            f"Falha ao ler experimental '{caminho_exp}': {e}"
        )

    cfast_embutido = serie.pop('cfast', None)
    resultado['experimental'] = {
        'label': exp.get('label', 'Experimental'),
        'cor':   exp.get('cor', 'magenta'),
        'serie': serie,
    }
    if cfast_embutido is not None and val.get('usar_cfast_embutido', False):
        resultado['cfast'].append({
            'label': cfast_embutido['fonte'],
            'cor':   'black',
            'serie': cfast_embutido,
        })

    # --- CFAST externos (opcional) ---
    for item in val.get('cfast', []) or []:
        caminho = _resolver_caminho(item['arquivo'])
        if not caminho.is_file():
            print(f"[AVISO] CFAST '{item.get('label', item['arquivo'])}' "
                  f"não encontrado em {caminho}. Pulando.")
            continue
        try:
            serie = ler_serie(caminho, formato=item.get('formato', 'auto'))
        except Exception as e:
            print(f"[AVISO] Falha ao ler '{caminho}': {e}. Pulando.")
            continue
        resultado['cfast'].append({
            'label': item.get('label', serie['fonte']),
            'cor':   item.get('cor', 'black'),
            'serie': serie,
        })

    # --- Modelos / saídas do simulador (opcional, mas exigimos pelo menos um) ---
    for item in val.get('modelos', []) or []:
        caminho = _resolver_caminho(item['arquivo'])
        if not caminho.is_file():
            print(f"[AVISO] Modelo '{item.get('label', item['arquivo'])}' "
                  f"não encontrado em {caminho}. Pulando.")
            continue
        try:
            serie = ler_serie(caminho, formato='simulador')
        except Exception as e:
            print(f"[AVISO] Falha ao ler '{caminho}': {e}. Pulando.")
            continue
        resultado['modelos'].append({
            'label': item.get('label', 'Modelo'),
            'cor':   item.get('cor', 'red'),
            'serie': serie,
        })

    if not resultado['modelos']:
        raise ErroCenario(
            "Nenhum modelo do simulador foi carregado. Gere ao menos um CSV "
            "com `python main.py --cenario ... --saida csv` antes de validar."
        )

    return resultado

# ---------------------------------------------------------------------------
# Plotagem
# ---------------------------------------------------------------------------
def plotar_painel(cfg, dados, t_grid, t_max, saida_png):
    nome_cenario = cfg.get('nome', 'cenário')
    referencia   = cfg.get('referencia', '').strip()
    descricao    = cfg.get('descricao', '').strip()

    titulo = f"Validação — {nome_cenario}"
    if referencia:
        titulo += f"\n({referencia})"

    fig, axs = plt.subplots(5, 1, figsize=(11, 14), sharex=True)

    def _plotar(ax, chave, ylabel, escala=1.0, titulo_ax=''):
        # Experimental
        if dados['experimental']:
            s = dados['experimental']
            y = interpolar_em_grid(s['serie'], t_grid).get(chave)
            if y is not None:
                ax.plot(t_grid, y * escala, '.', color=s['cor'],
                        markersize=2, label=s['label'], zorder=4)
        # CFAST
        for s in dados['cfast']:
            y = interpolar_em_grid(s['serie'], t_grid).get(chave)
            if y is not None:
                ax.plot(t_grid, y * escala, 'o', color=s['cor'],
                        markersize=2, markerfacecolor='none',
                        label=s['label'], zorder=3)
        # Modelos
        for s in dados['modelos']:
            y = interpolar_em_grid(s['serie'], t_grid).get(chave)
            if y is not None:
                ax.plot(t_grid, y * escala, '-', color=s['cor'],
                        linewidth=1.5, label=s['label'], zorder=2)

        ax.set_ylabel(ylabel)
        ax.set_title(titulo_ax, loc='left', fontsize=10)
        ax.grid(True, alpha=0.3)
        ax.legend(fontsize=8, loc='best')

    _plotar(axs[0], 'HRR_kW',  'HRR (kW)',       titulo_ax='Taxa de liberação de calor')
    _plotar(axs[1], 'T_up_C',  'T (°C)',         titulo_ax='Temperatura da camada superior')
    _plotar(axs[2], 'Z_int_m', 'Z_int (m)',      titulo_ax='Altura da interface')
    _plotar(axs[3], 'O2_vol',  'O₂ (fração vol.)', escala=100.0,
            titulo_ax='Oxigênio (base volumétrica)')
    _plotar(axs[4], 'CO2_vol', 'CO₂ (fração vol.)', escala=100.0,
            titulo_ax='Dióxido de carbono (base volumétrica)')

    axs[4].set_xlabel('Tempo (s)')
    axs[0].set_xlim(0, t_max)

    fig.suptitle(titulo, fontweight='bold', fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.97])
    fig.savefig(saida_png, dpi=200)
    print(f"Figura salva em: {saida_png}")
    plt.close(fig)


# ---------------------------------------------------------------------------
# Comparação numérica
# ---------------------------------------------------------------------------
def _estatisticas_referencia(exp, t_grid):
    """Imprime um resumo das séries experimentais no grid de comparação."""
    print("\n=== REFERÊNCIA EXPERIMENTAL ===")
    grid_exp = interpolar_em_grid(exp['serie'], t_grid)
    print(f"{'Série':<15} {'n':>5} {'min':>10} {'max':>10} {'média':>10} {'desv.p.':>10}")
    print('-' * 65)
    for chave, unidade in [('T_up_C', '°C'),
                            ('Z_int_m', 'm'),
                            ('O2_vol',  'frac vol'),
                            ('CO2_vol', 'frac vol'),
                            ('HRR_kW',  'kW')]:
        y = grid_exp.get(chave)
        if y is None:
            continue
        m = np.isfinite(y)
        if m.sum() == 0:
            continue
        print(f"{chave:<15} {int(m.sum()):>5d} "
              f"{np.nanmin(y):>10.4f} {np.nanmax(y):>10.4f} "
              f"{np.nanmean(y):>10.4f} {np.nanstd(y):>10.4f}   [{unidade}]")
    print()


def comparar(dados, t_grid):
    exp = dados['experimental']
    if exp is None:
        print("[AVISO] Sem dados experimentais — pulando métricas.")
        return

    # Bloco descritivo da referência
    _estatisticas_referencia(exp, t_grid)

    grid_exp = interpolar_em_grid(exp['serie'], t_grid)

    print("\n=== MÉTRICAS vs EXPERIMENTAL ===\n")
    for chave, unidade in [('T_up_C', '°C'),
                            ('Z_int_m', 'm'),
                            ('O2_vol',  'fração vol.'),
                            ('CO2_vol', 'fração vol.')]:
        tabela = []
        ref = grid_exp.get(chave)
        if ref is None:
            continue
        for s in dados['modelos'] + dados['cfast']:
            y = interpolar_em_grid(s['serie'], t_grid).get(chave)
            if y is None:
                continue
            m = metricas(y, ref, f"{s['label']} — {chave} [{unidade}]")
            if m:
                tabela.append(m)
        if tabela:
            imprimir_metricas(tabela)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def main():
    parser = argparse.ArgumentParser(description="Validação do simulador")
    parser.add_argument('--cenario', '-c', required=True,
                        help='Nome do cenário (sem extensão) em ./cenarios/')
    parser.add_argument('--t-max', type=float, default=None,
                        help='Tempo máximo nos gráficos (default: t_final do cenário)')
    parser.add_argument('--t-step', type=float, default=1.0,
                        help='Passo do grid de interpolação [s]')
    parser.add_argument('--saida', default=None,
                        help='Nome do PNG de saída (default: validacao_<cenario>.png)')
    args = parser.parse_args()

    try:
        caminho = localizar_cenario(args.cenario)
        cfg = carregar_cenario(caminho)
    except ErroCenario as e:
        print(f"[ERRO] {e}", file=sys.stderr)
        sys.exit(1)

    t_final = cfg['simulacao']['t_final']
    t_max = args.t_max if args.t_max is not None else t_final
    t_grid = np.arange(0.0, t_max + args.t_step, args.t_step)

    print(f"Cenário: {cfg['nome']}")
    print(f"  Referência: {cfg.get('referencia', '—')}")
    print(f"  t_final (YAML): {t_final} s   |   t_max (plot): {t_max} s")

    dados = carregar_series(cfg)

    n_exp = 1 if dados['experimental'] else 0
    print(f"  Fontes: experimental={n_exp}, CFAST={len(dados['cfast'])}, "
          f"modelos={len(dados['modelos'])}")

    saida = args.saida or f"validacao_{cfg['nome']}.png"
    plotar_painel(cfg, dados, t_grid, t_max, saida)
    comparar(dados, t_grid)


if __name__ == "__main__":
    main()
