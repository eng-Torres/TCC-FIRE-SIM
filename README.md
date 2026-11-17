# Simulador de Incêndios em Compartimentos (Duas Zonas)

## 1. Descrição

Simulador computacional de incêndios em compartimentos baseado no modelo termofluidodinâmico de duas zonas em Python, com acoplamento térmico de paredes em 1D e ventilação natural por plano neutro. O código oferece suporte a três correlações de pluma térmica (`heskestad`, `mccaffrey` e `zukoski`), combustão de múltiplos objetos com ignição sequencial e rotina automatizada de validação contra dados de referência.

## 2. Requisitos e Instalação

Requer Python 3.10+ e as bibliotecas científicas abaixo:

```bash
pip install numpy scipy pandas matplotlib pyyaml
```

*(Alternativamente, `pip install -r requirements.txt`.)*

## 3. Como Rodar o Simulador (`main.py`)

Os cenários ficam declarados em `cenarios/` (`quarto_padrao.yaml` e `teste_babrauskas.yaml`). Para listar opções ou executar uma simulação escolhendo a correlação de pluma e direcionando o CSV de saída para a pasta `resultados/`:

```bash
# Listar cenários e materiais disponíveis
python main.py --listar
python main.py --listar-materiais

# Rodar o cenário Quarto Padrão (cenarios/quarto_padrao.yaml)
python main.py --cenario quarto_padrao --pluma heskestad --saida csv --nome-saida resultados/resultados_quarto_padrao_heskestad.csv

# Rodar o cenário Teste 1 de Babrauskas (cenarios/teste_babrauskas.yaml)
python main.py --cenario teste_babrauskas --pluma heskestad --saida csv --nome-saida resultados/resultados_teste1_Babrauskas_heskestad.csv
```

Use `--saida ambos` para exportar o `.csv` e exibir na tela os gráficos transientes de camada superior, interface de fumaça, temperatura e espécies (`--pluma` aceita `heskestad`, `mccaffrey` ou `zukoski`). Ao usar `--nome-saida resultados/...`, as séries temporais `.csv` são salvas na pasta `resultados/`, e os painéis comparativos `.png` da etapa seguinte ficam conforme configurado no YAML (ou redirecionados via `--saida`).

## 4. Como Rodar a Validação (`valid_cfast.py`)

Com os arquivos CSV das três plumas gerados na pasta `resultados/`, execute a validação inter-modelos (contra CFAST 7.7) e a validação experimental (contra `1Room1.csv` do NIST):

```bash
python valid_cfast.py --config config_validacao_quarto.yaml
python valid_cfast.py --config config_validacao_babrauskas.yaml
```

O primeiro comando compara as três plumas contra o CFAST 7.7 e gera o gráfico comparativo `fig_quarto_padrao_cfast.png` e a tabela de métricas estatísticas `metricas_quarto_padrao_cfast.csv` (RMSE, MAE, erro máximo, viés, Pearson *r* e ΔExt). O segundo confronta as simulações com os dados experimentais do NIST (`1Room1.csv`), salvando o painel `fig_babrauskas_valid_cfast.png` e as métricas em `metricas_babrauskas_valid_cfast.csv` (use `--saida resultados/<arquivo>.png` caso deseje gravar o gráfico diretamente em `resultados/`).
