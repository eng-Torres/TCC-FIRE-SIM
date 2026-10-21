#!/usr/bin/env bash
# rodar_plumas.sh
for pluma in mccaffrey heskestad zukoski; do
    echo "=== Rodando $pluma ==="
    python main.py --cenario teste_babrauskas --pluma "$pluma"  --nome-saida cenarios/resultados_teste1_Babrauskas_"$pluma".csv
done
python valid.py --cenario teste_babrauskas --saida fig_teste1.png