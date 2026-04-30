"""
Módulo de Geometria do Simulador de Incêndio (Modelo de Duas Zonas)
Responsável por armazenar as dimensões do compartimento, objetos e aberturas (vents).
"""

def carregar_geometria_padrao():
    """
    Retorna a geometria padrão de calibração: um quarto de 4x4x2.8m,
    contendo uma cama, um guarda-roupa e uma porta aberta.
    """
    
    # 1. O Mestre: A Sala
    # Guardamos as dimensões e já deixamos a área e o volume pré-calculados
    # para poupar processamento nas equações da física depois.
    sala = {
        "comprimento_x": 4.0,
        "profundidade_y": 4.0,
        "altura_z": 2.8,
        "area_piso": 4.0 * 4.0,
        "volume_total": 4.0 * 4.0 * 2.8
    }

    # 2. A Lista de Inventário: Móveis (Sólidos Combustíveis)
    # Cada dicionário é um móvel com sua posição (ponto inferior esquerdo) e tamanho.
    objetos = [
        {
            "nome": "Guarda-Roupa",
            "posicao": {"x": 0.0, "y": 1.0, "z": 0.0},
            "dimensoes": {"x": 1.5, "y": 0.6, "z": 2.2},
            "tipo_material": "madeira" # Gancho para a física ler depois
        },
        {
            "nome": "Cama de Casal",
            "posicao": {"x": 2.6, "y": 2.0, "z": 0.0},
            "dimensoes": {"x": 1.4, "y": 2.0, "z": 0.5},
            "tipo_material": "espuma_poliuretano"
        }
    ]

    # 3. As Aberturas: Respiros de Oxigênio (Vents)
    # A profundidade (y) é 0 porque a porta é apenas uma fenda bidimensional na parede.
    aberturas = [
        {
            "nome": "Porta Principal",
            "posicao": {"x": 1.6, "y": 0.0, "z": 0.0},
            "dimensoes": {"x": 0.8, "y": 0.0, "z": 2.1} 
        }
    ]

    return sala, objetos, aberturas