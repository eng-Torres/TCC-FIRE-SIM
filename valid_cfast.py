"""
Script genérico de validação: Simulador Python vs. CFAST (NIST) orientado por YAML.

Uso:
    python valid_cfast.py
    python valid_cfast.py --config config_validacao.yaml
    python valid_cfast.py --config config_validacao.yaml --saida comparacao_cfast.png

Funcionalidades:
- Leitura guiada por arquivo de configuração YAML (sem caminhos ou colunas fixas no código).
- Autodetecção dos cabeçalhos de 4 linhas do CFAST 7.x (incluindo linha de unidades).
- Conversão automática de 'mol %' do CFAST para fração molar/volumétrica [0–1] (ou via `cfast_scale`).
- Mapeamento inteligente de apelidos de colunas do simulador (`T_u` -> `T_u(C)`,
  `Z_int` -> `Z_int(m)`, `T_teto`/`T_w` -> `T_w(C)`, `O2_vol`, `CO2_vol`, `CO_vol`).
- Cálculo de métricas estatísticas (RMSE, MAE, Erro Máximo, Viés, Erro Relativo no Pico/Vale
  segundo ASTM E1355 e Correlação de Pearson r).
- Exportação de painel gráfico de alta resolução e tabela de métricas em CSV.
"""

import argparse
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import yaml

# Massas molares [g/mol] para conversão fração mássica (Y_i) -> fração volumétrica (X_i)
M_AR = 28.97
M_O2 = 32.00
M_CO2 = 44.01
M_CO = 28.01

# Apelidos de colunas do CSV do simulador Python (outputs.py)
ALIASES_SIMULADOR = {
    "t": "Tempo(s)",
    "time": "Tempo(s)",
    "tempo": "Tempo(s)",
    "z_int": "Z_int(m)",
    "hgt": "Z_int(m)",
    "t_u": "T_u(C)",
    "ult": "T_u(C)",
    "t_w": "T_w(C)",
    "t_teto": "T_w(C)",
    "t_parede": "T_w(C)",
    "ceilt": "T_w(C)",
    "m_dot_in": "m_dot_in(kg/s)",
    "m_dot_out": "m_dot_out(kg/s)",
}


# ---------------------------------------------------------------------------
# Leitura e normalização do CSV do Simulador Python
# ---------------------------------------------------------------------------
def carregar_csv_simulador(caminho: Path) -> pd.DataFrame:
    """
    Lê o CSV exportado por src/outputs.py e adiciona colunas derivadas
    (O2_vol, CO2_vol, CO_vol, O2_pct, CO2_pct) e apelidos padronizados.
    """
    if not caminho.is_file():
        raise FileNotFoundError(f"Arquivo do simulador não encontrado: {caminho}")

    df = pd.read_csv(caminho, comment="#")
    df.columns = [c.strip() for c in df.columns]

    # Conversão de fração mássica para fração volumétrica (mol/mol)
    if "Y_O2" in df.columns and "O2_vol" not in df.columns:
        df["O2_vol"] = np.clip(df["Y_O2"].to_numpy(float), 0.0, 1.0) * (M_AR / M_O2)
        df["O2_pct"] = df["O2_vol"] * 100.0
    if "Y_CO2" in df.columns and "CO2_vol" not in df.columns:
        df["CO2_vol"] = np.clip(df["Y_CO2"].to_numpy(float), 0.0, 1.0) * (M_AR / M_CO2)
        df["CO2_pct"] = df["CO2_vol"] * 100.0
    if "Y_CO" in df.columns and "CO_vol" not in df.columns:
        df["CO_vol"] = np.clip(df["Y_CO"].to_numpy(float), 0.0, 1.0) * (M_AR / M_CO)
        df["CO_pct"] = df["CO_vol"] * 100.0

    return df


def extrair_coluna_simulador(df: pd.DataFrame, nome_col: str) -> np.ndarray:
    """Extrai uma coluna do DataFrame do simulador aceitando nomes exatos ou apelidos."""
    if nome_col in df.columns:
        return pd.to_numeric(df[nome_col], errors="coerce").to_numpy(float)

    chave = nome_col.strip().lower()
    alvo = ALIASES_SIMULADOR.get(chave)
    if alvo and alvo in df.columns:
        return pd.to_numeric(df[alvo], errors="coerce").to_numpy(float)

    # Busca case-insensitive ou por prefixo antes de '('
    for col in df.columns:
        prefixo = col.split("(")[0].strip().lower()
        if col.lower() == chave or prefixo == chave:
            return pd.to_numeric(df[col], errors="coerce").to_numpy(float)

    raise KeyError(
        f"Coluna '{nome_col}' não encontrada no CSV do simulador. "
        f"Colunas disponíveis: {list(df.columns)}"
    )


# ---------------------------------------------------------------------------
# Leitura dos CSVs do CFAST 7.x ou Experimentais Indexados (ex.: 1Room1.csv)
# ---------------------------------------------------------------------------
def carregar_csv_cfast(
    caminho: Path,
    skiprows: int | list[int] | None = None,
    header_none: bool = False,
) -> tuple[pd.DataFrame, dict[int | str, str]]:
    """
    Lê um arquivo CSV de referência:
    1) Modo indexado/legado (`header_none=True`, ex.: `1Room1.csv` do NIST com `skiprows=2, header=None`):
       Usa `pd.read_csv(caminho, skiprows=skiprows, header=None)` e acessa colunas por índice inteiro.
    2) Modo padrão CFAST 7.x (compartments, walls, vents, zone, etc.):
       Possui 4 linhas de cabeçalho (ID curto, descrição, alvo, unidade).
    Retorna (df_numerico, dict_unidades_por_coluna).
    """
    if not caminho.is_file():
        raise FileNotFoundError(f"Arquivo de referência/CFAST não encontrado: {caminho}")

    if header_none:
        sr = 2 if skiprows is None else skiprows
        df = pd.read_csv(caminho, skiprows=sr, header=None)
        unidades: dict[int | str, str] = {}
        # Em 1Room1.csv (skiprows=2), a linha 1 do df (4ª linha do arquivo) contém as unidades
        if len(df) >= 2:
            for col_idx in df.columns:
                val_u = str(df.iloc[1, col_idx]).strip() if pd.notna(df.iloc[1, col_idx]) else ""
                unidades[int(col_idx)] = val_u
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors="coerce")
        return df, unidades

    with open(caminho, encoding="utf-8", errors="ignore") as f:
        primeiras = [f.readline() for _ in range(5)]

    sep = ";" if primeiras[0].count(";") > primeiras[0].count(",") else ","

    cabecalho_ids = [c.strip() for c in primeiras[0].strip().split(sep)]
    unidades_raw = [u.strip() for u in primeiras[3].strip().split(sep)] if len(primeiras) >= 4 else []

    unidades = {}
    for i, col_id in enumerate(cabecalho_ids):
        unidades[col_id] = unidades_raw[i] if i < len(unidades_raw) else ""

    df = pd.read_csv(caminho, sep=sep, header=0, skiprows=[1, 2, 3])
    df.columns = [c.strip() for c in df.columns]
    for col in df.columns:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    return df, unidades


def extrair_coluna_cfast(
    df: pd.DataFrame,
    unidades: dict[int | str, str],
    nome_col: str | int,
    sim_col: str = "",
    escala_manual: float | None = None,
    offset_manual: float | None = None,
) -> tuple[np.ndarray, str]:
    """
    Extrai a coluna `nome_col` (índice inteiro para CSV sem cabeçalho como 1Room1.csv,
    ou nome/soma de colunas separadas por '+' para CFAST 7.x) e aplica conversão de unidades.
    """
    if isinstance(nome_col, int) or (isinstance(nome_col, str) and nome_col.strip().isdigit()):
        idx = int(nome_col)
        if idx not in df.columns:
            raise KeyError(
                f"Índice de coluna {idx} não encontrado no CSV ({df.shape[1]} colunas disponíveis)."
            )
        valores = df[idx].to_numpy(float).copy()
        unidade_detectada = str(unidades.get(idx, ""))
    else:
        partes = [p.strip() for p in str(nome_col).split("+")]
        valores = np.zeros(len(df), dtype=float)
        unidade_detectada = ""

        for p in partes:
            if p not in df.columns:
                raise KeyError(
                    f"Coluna '{p}' não encontrada no CSV do CFAST. "
                    f"Colunas disponíveis: {list(df.columns)}"
                )
            valores += df[p].to_numpy(float)
            if not unidade_detectada:
                unidade_detectada = str(unidades.get(p, ""))

    # Conversão automática de temperatura em Kelvin para Celsius (ex.: TEMPUP no 1Room1.csv)
    sim_col_norm = (sim_col or "").strip().lower()
    if offset_manual is not None:
        valores = valores + float(offset_manual)
    elif (
        unidade_detectada.strip().upper() == "K"
        or (sim_col_norm in ("t_u", "t_u(c)", "ult") and np.nanmin(valores) > 200.0)
    ):
        valores = valores - 273.15

    if escala_manual is not None:
        valores = valores * float(escala_manual)
    else:
        # Se o arquivo exporta em '%' ou 'mol %' (0-100) e o simulador compara em fração ('*_vol')
        eh_mol_pct = "mol %" in unidade_detectada.lower() or "%" in unidade_detectada.lower()
        quer_fracao = sim_col_norm.endswith("_vol")
        if (eh_mol_pct and quer_fracao) or (quer_fracao and np.nanmax(valores) > 1.0):
            valores = valores / 100.0

    return valores, unidade_detectada


# ---------------------------------------------------------------------------
# Interpolação e Métricas Estatísticas (ASTM E1355)
# ---------------------------------------------------------------------------
def interpolar_serie(t_orig: np.ndarray, y_orig: np.ndarray, t_grid: np.ndarray) -> np.ndarray:
    """Interpola linearmente (t_orig, y_orig) na malha `t_grid`, ignorando NaNs."""
    mask = np.isfinite(t_orig) & np.isfinite(y_orig)
    if mask.sum() < 2:
        return np.full_like(t_grid, np.nan, dtype=float)
    return np.interp(t_grid, t_orig[mask], y_orig[mask], left=np.nan, right=np.nan)


def calcular_metricas(
    y_sim: np.ndarray,
    y_ref: np.ndarray,
    nome_var: str,
    modelo: str,
    escala_extremo: float | None = None,
) -> dict:
    """
    Calcula RMSE, MAE, Erro Máximo, Viés, Correlação de Pearson (r) e
    Erro Relativo de Pico/Vale (%) conforme diretrizes da ASTM E1355.

    Quando `escala_extremo` é fornecida (ex.: H_sala para Z_int), o erro
    relativo de extremo é normalizado por essa escala característica em vez
    de |ext_ref|, evitando inflação numérica quando ext_ref -> 0.
    """
    mask = np.isfinite(y_sim) & np.isfinite(y_ref)
    if mask.sum() == 0:
        return {}

    ys = y_sim[mask]
    yr = y_ref[mask]
    diff = ys - yr

    rmse = float(np.sqrt(np.mean(diff ** 2)))
    mae = float(np.mean(np.abs(diff)))
    max_err = float(np.max(np.abs(diff)))
    bias = float(np.mean(diff))

    # Correlação de Pearson
    if np.std(ys) > 1e-12 and np.std(yr) > 1e-12:
        r_pearson = float(np.corrcoef(ys, yr)[0, 1])
    else:
        r_pearson = np.nan

    # Extremo característico (vale para altura/O2, pico para temperatura/CO2/HRR)
    # Detecta se a variável diminui em relação ao valor inicial (ex: Z_int, O2)
    cai_no_incendio = (yr[0] - np.min(yr)) > (np.max(yr) - yr[0])
    if cai_no_incendio:
        ext_sim = float(np.min(ys))
        ext_ref = float(np.min(yr))
        tipo_ext = "mín"
    else:
        ext_sim = float(np.max(ys))
        ext_ref = float(np.max(yr))
        tipo_ext = "máx"

    denom_ext = float(escala_extremo) if (escala_extremo is not None and escala_extremo > 0.0) else abs(ext_ref)
    erro_rel_ext_pct = (
        float((ext_sim - ext_ref) / denom_ext * 100.0)
        if denom_ext > 1e-12
        else np.nan
    )

    return {
        "variavel": nome_var,
        "modelo": modelo,
        "n": int(mask.sum()),
        "rmse": rmse,
        "mae": mae,
        "max_err": max_err,
        "bias": bias,
        "r_pearson": r_pearson,
        "tipo_extremo": tipo_ext,
        "extremo_sim": ext_sim,
        "extremo_cfast": ext_ref,
        "erro_rel_extremo_pct": erro_rel_ext_pct,
        "escala_extremo": float(escala_extremo) if escala_extremo is not None else np.nan,
    }


def imprimir_relatorio(lista_metricas: list[dict], cenario: str):
    """Imprime a tabela de métricas no terminal."""
    print(f"\n{'=' * 108}")
    print(f" RELATÓRIO DE VALIDAÇÃO CONTRA O CFAST — CENÁRIO: {cenario.upper()}")
    print(f"{'=' * 108}")
    cab = (
        f"{'Variável':<24} {'Modelo':<22} {'n':>4} "
        f"{'RMSE':>9} {'MAE':>9} {'Máx':>9} {'Viés':>9} "
        f"{'r':>7} {'Ext.Sim':>9} {'Ext.CFAST':>9} {'ΔExt(%)':>8}"
    )
    print(cab)
    print("-" * 108)
    for m in lista_metricas:
        print(
            f"{m['variavel']:<24} {m['modelo']:<22} {m['n']:>4d} "
            f"{m['rmse']:>9.4f} {m['mae']:>9.4f} {m['max_err']:>9.4f} {m['bias']:>+9.4f} "
            f"{m['r_pearson']:>7.4f} {m['extremo_sim']:>9.4f} {m['extremo_cfast']:>9.4f} "
            f"{m['erro_rel_extremo_pct']:>+7.2f}%"
        )
    print("=" * 108 + "\n")


# ---------------------------------------------------------------------------
# Execução Principal e Plotagem
# ---------------------------------------------------------------------------
def executar_validacao(caminho_config: Path, saida_png_cli: str | None = None):
    if not caminho_config.is_file():
        raise FileNotFoundError(f"Arquivo de configuração não encontrado: {caminho_config}")

    base_dir = caminho_config.resolve().parent
    with open(caminho_config, encoding="utf-8") as f:
        cfg = yaml.safe_load(f)

    cenario = cfg.get("cenario", "cenario")
    altura_sala_cfg = cfg.get("altura_sala")
    titulo_geral = cfg.get("titulo", f"Validação Simulador vs. CFAST — {cenario}")
    pasta_cfast = base_dir / cfg["caminho_cfast"]
    tempo_max = float(cfg.get("tempo_max", 600.0))
    dt_grid = float(cfg.get("dt_grid", 1.0))
    t_grid = np.arange(0.0, tempo_max + 0.5 * dt_grid, dt_grid)

    # Lista de modelos do simulador (suporta 1 arquivo ou múltiplos)
    modelos_cfg = cfg.get("modelos_simulador")
    if not modelos_cfg:
        caminho_sim = base_dir / cfg["caminho_simulador"]
        modelos_cfg = [{
            "arquivo": str(caminho_sim),
            "label": "Simulador Python",
            "cor": "crimson",
        }]

    dados_modelos = []
    for item in modelos_cfg:
        p_sim = Path(item["arquivo"])
        if not p_sim.is_absolute():
            p_sim = base_dir / p_sim
        df_sim = carregar_csv_simulador(p_sim)
        t_sim = extrair_coluna_simulador(df_sim, "Tempo(s)")
        dados_modelos.append({
            "label": item.get("label", p_sim.stem),
            "cor": item.get("cor", "crimson"),
            "df": df_sim,
            "t": t_sim,
        })

    variaveis = cfg.get("variaveis", {})
    if not variaveis:
        raise ValueError("Nenhuma variável definida na seção 'variaveis' do YAML.")

    # Cache de arquivos de referência/CFAST já lidos
    cache_cfast: dict[tuple[Path, bool, int | None], tuple[pd.DataFrame, dict[int | str, str]]] = {}

    n_vars = len(variaveis)
    plt.style.use("seaborn-v0_8-whitegrid" if "seaborn-v0_8-whitegrid" in plt.style.available else "default")
    fig, axs = plt.subplots(n_vars, 1, figsize=(11, 2.9 * n_vars), sharex=True)
    if n_vars == 1:
        axs = [axs]

    todas_metricas = []

    for ax, (chave_var, vcfg) in zip(axs, variaveis.items()):
        # Resolve arquivo de referência: aceita pasta_cfast como diretório ou arquivo direto
        cfast_file_cfg = vcfg.get("cfast_file")
        if pasta_cfast.is_file():
            arq_cfast = pasta_cfast
        elif cfast_file_cfg:
            cand = pasta_cfast / cfast_file_cfg
            if cand.is_file():
                arq_cfast = cand
            elif (base_dir / cfast_file_cfg).is_file():
                arq_cfast = base_dir / cfast_file_cfg
            elif (base_dir / "cenarios" / cfast_file_cfg).is_file():
                arq_cfast = base_dir / "cenarios" / cfast_file_cfg
            else:
                arq_cfast = cand
        else:
            cand_cen = base_dir / "cenarios" / pasta_cfast.name
            arq_cfast = cand_cen if cand_cen.is_file() else pasta_cfast

        col_ref = vcfg["cfast_col"]
        col_eh_indice = isinstance(col_ref, int) or (
            isinstance(col_ref, str) and col_ref.strip().isdigit()
        )
        header_none = bool(
            vcfg.get("cfast_header_none", cfg.get("cfast_header_none", col_eh_indice))
        )
        skiprows_cfg = vcfg.get("cfast_skiprows", cfg.get("cfast_skiprows", 2 if header_none else None))

        chave_cache = (arq_cfast, header_none, skiprows_cfg if isinstance(skiprows_cfg, int) else None)
        if chave_cache not in cache_cfast:
            cache_cfast[chave_cache] = carregar_csv_cfast(
                arq_cfast, skiprows=skiprows_cfg, header_none=header_none
            )
        df_cfast, unids_cfast = cache_cfast[chave_cache]

        time_col = vcfg.get("cfast_time_col", cfg.get("cfast_time_col", 86 if header_none else "Time"))
        if isinstance(time_col, int) or (isinstance(time_col, str) and time_col.strip().isdigit()):
            t_cfast = df_cfast[int(time_col)].to_numpy(float)
        else:
            t_cfast = df_cfast[str(time_col)].to_numpy(float)

        sim_col_cfg = vcfg.get("sim_col") or ""
        y_cfast_raw, unid_str = extrair_coluna_cfast(
            df_cfast,
            unids_cfast,
            col_ref,
            sim_col=sim_col_cfg,
            escala_manual=vcfg.get("cfast_scale"),
            offset_manual=vcfg.get("cfast_offset"),
        )
        y_cfast_grid = interpolar_serie(t_cfast, y_cfast_raw, t_grid)

        label_ref = vcfg.get(
            "cfast_label",
            cfg.get("cfast_label", f"Ref/CFAST (col {col_ref})" if header_none else f"CFAST ({col_ref})"),
        )
        # Plota curva de referência (CFAST ou Experimental 1Room1.csv)
        ax.plot(
            t_grid,
            y_cfast_grid,
            "k--",
            linewidth=2.0,
            label=label_ref,
            zorder=4,
        )

        # Plota cada modelo do simulador e calcula métricas (se a coluna existir no simulador)
        resumo_box = []
        if sim_col_cfg:
            for mod in dados_modelos:
                try:
                    y_sim_raw = extrair_coluna_simulador(mod["df"], sim_col_cfg)
                except KeyError:
                    continue
                if "sim_scale" in vcfg and vcfg["sim_scale"] is not None:
                    y_sim_raw = y_sim_raw * float(vcfg["sim_scale"])
                y_sim_grid = interpolar_serie(mod["t"], y_sim_raw, t_grid)

                escala_ext_var = vcfg.get("escala_extremo")
                if (
                    escala_ext_var is None
                    and altura_sala_cfg is not None
                    and (sim_col_cfg == "Z_int" or "interface" in chave_var.lower())
                ):
                    escala_ext_var = float(altura_sala_cfg)

                m = calcular_metricas(
                    y_sim_grid,
                    y_cfast_grid,
                    chave_var,
                    mod["label"],
                    escala_extremo=escala_ext_var,
                )
                if m:
                    todas_metricas.append(m)
                    resumo_box.append(
                        f"{mod['label']}: RMSE={m['rmse']:.3g} | MAE={m['mae']:.3g} | Viés={m['bias']:+.3g}"
                    )

                ax.plot(
                    t_grid,
                    y_sim_grid,
                    "-",
                    color=mod["cor"],
                    linewidth=1.8,
                    label=f"{mod['label']} ({sim_col_cfg})",
                    zorder=3,
                )

        titulo_sub = vcfg.get("titulo", chave_var.replace("_", " ").title())
        ax.set_title(titulo_sub, loc="left", fontsize=10.5, fontweight="bold")
        ax.set_ylabel(vcfg.get("ylabel", chave_var))
        ax.grid(True, alpha=0.35)
        ax.legend(loc="upper right", fontsize=8.5, frameon=True)

        if resumo_box:
            ax.text(
                0.015,
                0.92,
                "\n".join(resumo_box),
                transform=ax.transAxes,
                fontsize=8.0,
                verticalalignment="top",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.85, edgecolor="silver"),
            )

    axs[-1].set_xlabel("Tempo [s]", fontsize=10.5)
    axs[0].set_xlim(0.0, tempo_max)

    fig.suptitle(titulo_geral, fontsize=12.5, fontweight="bold")
    fig.tight_layout(rect=[0, 0, 1, 0.97])

    saida_png = saida_png_cli or cfg.get("saida_png", f"validacao_cfast_{cenario}.png")
    caminho_png = base_dir / saida_png if not Path(saida_png).is_absolute() else Path(saida_png)
    fig.savefig(caminho_png, dpi=200)
    plt.close(fig)

    imprimir_relatorio(todas_metricas, cenario)
    print(f"Gráfico comparativo salvo em: {caminho_png}")

    saida_csv = cfg.get("saida_metricas_csv")
    if saida_csv and todas_metricas:
        caminho_csv = base_dir / saida_csv if not Path(saida_csv).is_absolute() else Path(saida_csv)
        pd.DataFrame(todas_metricas).to_csv(caminho_csv, index=False, float_format="%.6f")
        print(f"Tabela de métricas exportada para: {caminho_csv}")


def main():
    parser = argparse.ArgumentParser(
        description="Validação genérica Simulador Python vs. CFAST via configuração YAML."
    )
    parser.add_argument(
        "--config",
        "-c",
        default="config_validacao.yaml",
        help="Caminho do arquivo YAML de configuração (default: config_validacao.yaml)",
    )
    parser.add_argument(
        "--saida",
        "-o",
        default=None,
        help="Override do caminho do arquivo PNG de saída",
    )
    args = parser.parse_args()

    try:
        executar_validacao(Path(args.config), saida_png_cli=args.saida)
    except Exception as e:
        print(f"[ERRO] {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
