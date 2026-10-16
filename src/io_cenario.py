"""
Leitura e validação de cenários (YAML ou JSON).

Suporta:
- Referência a materiais existentes em BANCO_MATERIAIS e MATERIAIS_CONSTRUCAO.
- Definição inline de materiais customizados no próprio arquivo de cenário.
- Compatibilidade retroativa com JSON (extensão .json).

"""

import json
from copy import deepcopy
from pathlib import Path

import yaml

from src.banco_materiais import BANCO_MATERIAIS
from src.banco_construcao import MATERIAIS_CONSTRUCAO

CENARIOS_DIR = Path(__file__).resolve().parent.parent / "cenarios"

PLUMAS_VALIDAS = ("heskestad", "zukoski", "mccaffrey")
MODOS_SAIDA = ("csv", "plot", "ambos")

CAMPOS_MAT_COMB = ("tempo", "hrr", "delta_Hc", "nu_o2", "y_co2", "y_co")
CAMPOS_MAT_CONSTR = ("k", "rho", "cp", "espessura", "h_conv", "emissividade")


class ErroCenario(Exception):
    """Erro de validação do cenário."""


# ---------------------------------------------------------------------------
# Leitura crua (YAML ou JSON)
# ---------------------------------------------------------------------------
def _ler_arquivo(caminho):
    """Lê um arquivo de cenário em YAML ou JSON, conforme a extensão."""
    sufixo = caminho.suffix.lower()

    if sufixo in (".yaml", ".yml"):
        with open(caminho, encoding="utf-8") as f:
            try:
                return yaml.safe_load(f)
            except yaml.YAMLError as e:
                raise ErroCenario(f"Erro ao parsear YAML: {e}")

    if sufixo == ".json":
        with open(caminho, encoding="utf-8") as f:
            try:
                return json.load(f)
            except json.JSONDecodeError as e:
                raise ErroCenario(f"Erro ao parsear JSON: {e}")

    raise ErroCenario(
        f"Extensão '{sufixo}' não suportada. Use .yaml, .yml ou .json"
    )


# ---------------------------------------------------------------------------
# Validação de materiais customizados
# ---------------------------------------------------------------------------
def _validar_material_combustivel(nome, mat):
    if not isinstance(mat, dict):
        raise ErroCenario(f"material combustível '{nome}' deve ser um dicionario")

    for k in CAMPOS_MAT_COMB:
        if k not in mat:
            raise ErroCenario(f"material combustível '{nome}': campo '{k}' ausente")

    tempo, hrr = mat["tempo"], mat["hrr"]
    if not isinstance(tempo, list) or not isinstance(hrr, list):
        raise ErroCenario(
            f"material combustível '{nome}': 'tempo' e 'hrr' devem ser listas"
        )
    if len(tempo) != len(hrr):
        raise ErroCenario(
            f"material combustível '{nome}': 'tempo' ({len(tempo)}) e "
            f"'hrr' ({len(hrr)}) têm tamanhos diferentes"
        )
    if len(tempo) < 2:
        raise ErroCenario(
            f"material combustível '{nome}': pelo menos 2 pontos são necessários"
        )
    for i in range(1, len(tempo)):
        if tempo[i] <= tempo[i - 1]:
            raise ErroCenario(
                f"material combustível '{nome}': 'tempo' deve ser estritamente "
                f"crescente (índice {i}: {tempo[i-1]} -> {tempo[i]})"
            )
    for i, h in enumerate(hrr):
        if h < 0:
            raise ErroCenario(
                f"material combustível '{nome}': 'hrr[{i}]' = {h} é negativo"
            )

    for k in ("delta_Hc", "nu_o2", "y_co2", "y_co"):
        if mat[k] <= 0:
            raise ErroCenario(f"material combustível '{nome}': '{k}' deve ser > 0")
    return mat


def _validar_material_construcao(nome, mat):
    if not isinstance(mat, dict):
        raise ErroCenario(f"material de construção '{nome}' deve ser um dicionario")

    for k in CAMPOS_MAT_CONSTR:
        if k not in mat:
            raise ErroCenario(f"material de construção '{nome}': campo '{k}' ausente")

    for k in ("k", "rho", "cp", "espessura", "h_conv"):
        if mat[k] <= 0:
            raise ErroCenario(f"material de construção '{nome}': '{k}' deve ser > 0")

    if not (0 < mat["emissividade"] <= 1):
        raise ErroCenario(
            f"material de construção '{nome}': 'emissividade' deve estar em (0, 1]"
        )
    return mat


# ---------------------------------------------------------------------------
# Validação estrutural
# ---------------------------------------------------------------------------
def _validar_sala(sala, mat_constr_validos):
    for k in ("comprimento_x", "profundidade_y", "altura_z"):
        if k not in sala:
            raise ErroCenario(f"sala: campo obrigatório '{k}' ausente")
        if not isinstance(sala[k], (int, float)) or sala[k] <= 0:
            raise ErroCenario(f"sala.{k} deve ser um número > 0")

    sala.setdefault("area_piso", sala["comprimento_x"] * sala["profundidade_y"])
    sala.setdefault("volume_total", sala["area_piso"] * sala["altura_z"])
    sala.setdefault("material_parede", "gesso_16mm")

    if sala["material_parede"] not in mat_constr_validos:
        raise ErroCenario(
            f"sala.material_parede '{sala['material_parede']}' não encontrado "
            f"no banco nem nos materiais customizados. "
            f"Disponíveis: {sorted(mat_constr_validos)}"
        )
    return sala


def _validar_aberturas(aberturas):
    if not aberturas:
        raise ErroCenario("Pelo menos uma abertura deve ser definida")
    for i, ab in enumerate(aberturas):
        for k in ("nome", "posicao", "dimensoes"):
            if k not in ab:
                raise ErroCenario(f"aberturas[{i}]: campo '{k}' ausente")
        for k in ("x", "y", "z"):
            if k not in ab["dimensoes"]:
                raise ErroCenario(f"aberturas[{i}].dimensoes.{k} ausente")
        if ab["dimensoes"]["x"] <= 0 or ab["dimensoes"]["z"] <= 0:
            raise ErroCenario(
                f"aberturas[{i}]: dimensoes.x e dimensoes.z devem ser > 0"
            )
    return aberturas


def _validar_objetos(objetos, mat_comb_validos):
    if not objetos:
        raise ErroCenario("Pelo menos um objeto deve ser definido")

    n_iniciais = 0
    for i, obj in enumerate(objetos):
        for k in ("nome_material", "area", "is_initial"):
            if k not in obj:
                raise ErroCenario(f"objetos[{i}]: campo '{k}' ausente")

        if obj["nome_material"] not in mat_comb_validos:
            raise ErroCenario(
                f"objetos[{i}].nome_material '{obj['nome_material']}' não "
                f"encontrado no banco nem nos materiais customizados. "
                f"Disponíveis: {sorted(mat_comb_validos)}"
            )
        if obj["area"] <= 0:
            raise ErroCenario(f"objetos[{i}].area deve ser > 0")

        if obj["is_initial"]:
            n_iniciais += 1
            obj.setdefault("T_ign", None)
        else:
            if "T_ign" not in obj or obj["T_ign"] is None:
                raise ErroCenario(
                    f"objetos[{i}]: objetos não-iniciais requerem T_ign"
                )

    if n_iniciais != 1:
        raise ErroCenario(
            f"Exatamente um objeto deve ter is_initial=true. "
            f"Encontrados: {n_iniciais}"
        )
    return objetos


def _validar_simulacao(sim):
    sim.setdefault("t_final", 600.0)
    sim.setdefault("modelo_pluma", "mccaffrey")
    sim.setdefault("max_step", 0.1)
    sim.setdefault("pontos_por_segundo", 10)
    sim.setdefault("n_nos_parede", 5)
    sim.setdefault("dt_print", 10.0)

    if sim["t_final"] <= 0:
        raise ErroCenario("simulacao.t_final deve ser > 0")
    if sim["modelo_pluma"] not in PLUMAS_VALIDAS:
        raise ErroCenario(
            f"simulacao.modelo_pluma '{sim['modelo_pluma']}' inválido. "
            f"Use uma de {PLUMAS_VALIDAS}."
        )
    if sim["max_step"] <= 0:
        raise ErroCenario("simulacao.max_step deve ser > 0")
    if sim["pontos_por_segundo"] < 1:
        raise ErroCenario("simulacao.pontos_por_segundo deve ser >= 1")
    if sim["n_nos_parede"] < 2:
        raise ErroCenario("simulacao.n_nos_parede deve ser >= 2")
    if sim["dt_print"] <= 0:
        raise ErroCenario("simulacao.dt_print deve ser > 0")
    return sim


def _validar_saida(saida):
    saida.setdefault("modo", "csv")
    saida.setdefault("nome_ficheiro", "resultados.csv")
    if saida["modo"] not in MODOS_SAIDA:
        raise ErroCenario(
            f"saida.modo '{saida['modo']}' inválido. Use um de {MODOS_SAIDA}."
        )
    return saida


# ---------------------------------------------------------------------------
# API pública
# ---------------------------------------------------------------------------
def carregar_cenario(caminho):
    """
    Carrega e valida um cenário. Retorna o dicionário `cfg` já normalizado.

    Estrutura retornada:
        nome, descricao, referencia,
        sala, aberturas, objetos,
        simulacao, saida,
        materiais_combustiveis, materiais_construcao
    """
    caminho = Path(caminho)
    if not caminho.is_file():
        raise ErroCenario(f"Arquivo de cenário não encontrado: {caminho}")

    cfg = _ler_arquivo(caminho)
    if not isinstance(cfg, dict):
        raise ErroCenario(
            f"Conteúdo do cenário deve ser um dicionario/objeto, "
            f"recebido {type(cfg).__name__}"
        )

    for k in ("sala", "aberturas", "objetos", "simulacao"):
        if k not in cfg:
            raise ErroCenario(f"Campo obrigatório '{k}' ausente")

    cfg.setdefault("nome", caminho.stem)
    cfg.setdefault("descricao", "")
    cfg.setdefault("referencia", "")
    cfg.setdefault("saida", {})
    cfg.setdefault("materiais_combustiveis", {})
    cfg.setdefault("materiais_construcao", {})

    for nome, mat in cfg["materiais_combustiveis"].items():
        _validar_material_combustivel(nome, mat)
    for nome, mat in cfg["materiais_construcao"].items():
        _validar_material_construcao(nome, mat)

    mat_comb_validos = set(BANCO_MATERIAIS) | set(cfg["materiais_combustiveis"])
    mat_constr_validos = set(MATERIAIS_CONSTRUCAO) | set(cfg["materiais_construcao"])

    cfg["sala"]      = _validar_sala(cfg["sala"], mat_constr_validos)
    cfg["aberturas"] = _validar_aberturas(cfg["aberturas"])
    cfg["objetos"]   = _validar_objetos(cfg["objetos"], mat_comb_validos)
    cfg["simulacao"] = _validar_simulacao(cfg["simulacao"])
    cfg["saida"]     = _validar_saida(cfg["saida"])
    return cfg


def registrar_materiais_customizados(cfg):
    """
    Copia materiais customizados do cenário para os bancos globais.

    Efeito colateral EXPLÍCITO: modifica BANCO_MATERIAIS e MATERIAIS_CONSTRUCAO.
    Chamar UMA vez por execução, antes de simular.

    Colisão: o material do cenário sobrescreve o do banco embutido, com aviso.
    """
    for nome, mat in cfg.get("materiais_combustiveis", {}).items():
        if nome in BANCO_MATERIAIS:
            print(f"[AVISO] material combustível '{nome}' sobrescreve o banco embutido.")
        BANCO_MATERIAIS[nome] = deepcopy(mat)

    for nome, mat in cfg.get("materiais_construcao", {}).items():
        if nome in MATERIAIS_CONSTRUCAO:
            print(f"[AVISO] material de construção '{nome}' sobrescreve o banco embutido.")
        MATERIAIS_CONSTRUCAO[nome] = deepcopy(mat)


def listar_cenarios():
    """Lista cenários disponíveis em ./cenarios/ (YAML ou JSON)."""
    if not CENARIOS_DIR.is_dir():
        return []
    extensoes = {".yaml", ".yml", ".json"}
    return sorted(
        p.stem for p in CENARIOS_DIR.iterdir()
        if p.is_file() and p.suffix.lower() in extensoes
    )


def listar_materiais():
    """Lista materiais do banco embutido. Útil para CLI."""
    return {
        "combustiveis": sorted(BANCO_MATERIAIS),
        "construcao": sorted(MATERIAIS_CONSTRUCAO),
    }


def localizar_cenario(nome):
    """
    Dado um nome sem extensão, retorna o Path do primeiro cenário
    encontrado com esse nome em ./cenarios/.
    """
    for ext in (".yaml", ".yml", ".json"):
        p = CENARIOS_DIR / f"{nome}{ext}"
        if p.is_file():
            return p
    raise ErroCenario(
        f"Cenário '{nome}' não encontrado em {CENARIOS_DIR}. "
        f"Use --listar para ver os disponíveis."
    )
