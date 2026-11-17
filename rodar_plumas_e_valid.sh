#!/usr/bin/env bash
# rodar_plumas.sh
for pluma in mccaffrey heskestad zukoski; do
    echo "=== Rodando $pluma ==="
    python main.py --cenario teste_babrauskas --pluma "$pluma"  --nome-saida cenarios/resultados_teste1_Babrauskas_"$pluma".csv
done
python valid_cfast.py --config config_validacao_babrauskas.yaml --saida fig_babrauskas_valid.png