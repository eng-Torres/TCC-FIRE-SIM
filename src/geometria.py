"""
Módulo de Geometria do Simulador de Incêndio (Modelo de Duas Zonas)
"""

def carregar_geometria_padrao():
    """Quarto 4×4×2.8 m com cama, guarda-roupa e porta aberta."""
    sala = {
        "comprimento_x": 4.0,
        "profundidade_y": 4.0,
        "altura_z": 2.8,
        "area_piso": 4.0 * 4.0,
        "volume_total": 4.0 * 4.0 * 2.8
    }

    objetos = [
        {
            "nome": "Guarda-Roupa",
            "posicao": {"x": 0.0, "y": 1.0, "z": 0.0},
            "dimensoes": {"x": 1.5, "y": 0.6, "z": 2.2},
            "tipo_material": "madeira",
            "T_ign": 573.15
        },
        {
            "nome": "Cama de Casal",
            "posicao": {"x": 2.6, "y": 2.0, "z": 0.0},
            "dimensoes": {"x": 1.4, "y": 2.0, "z": 0.5},
            "tipo_material": "espuma_poliuretano",
            "T_ign": 573.15
        },
        {
            "nome": "Tapete de Papelao",
            "posicao": {"x": 1.5, "y": 1.0, "z": 0.0},
            "dimensoes": {"x": 1.0, "y": 2.0, "z": 0.01},
            "tipo_material": "papelao",
            "T_ign": 373.15
        },
        {
            "nome": "cadeira",
            "posicao": {"x": 1.5, "y": 1.0, "z": 0.0},
            "dimensoes": {"x": 0.5, "y": 0.5, "z": 0.5},
            "tipo_material": "armchair_h",
            "T_ign": 473.15    



        }
    ]

    aberturas = [
        {
            "nome": "Porta Principal",
            "posicao": {"x": 1.6, "y": 0.0, "z": 0.0},
            "dimensoes": {"x": 0.8, "y": 0.0, "z": 2.1}
        }
    ]

    return sala, objetos, aberturas


def carregar_geometria_teste1_Babrauskas():
    """
    Geometria do Teste 1 de Babrauskas (1984).
    Sala: 2.26 × 3.94 × 2.31 m; janela 2.0 m (largura) × 1.13 m (altura).
    A√h = 2.43 m²·⁵ (Tabela 1).
    """
    sala = {
        "comprimento_x": 2.26,
        "profundidade_y": 3.94,
        "altura_z": 2.31,
        "area_piso": 2.26 * 3.94,
        "volume_total": 2.26 * 3.94 * 2.31
    }

    objetos = [
        {
            "nome": "Love Seat",
            "posicao": {"x": 1.13, "y": 1.0, "z": 0.0},
            "dimensoes": {"x": 1.5, "y": 0.8, "z": 0.9},
            "tipo_material": "teste1_hrr_medida",
            "T_ign": 573.15   # 300 °C — espuma de PU
        }
    ]

    aberturas = [
        {
            "nome": "Janela",
            "posicao": {"x": 0.13, "y": 0.0, "z": 0.0},   
            "dimensoes": {"x": 2.0, "y": 0.0, "z": 1.13}  
        }
    ]

    return sala, objetos, aberturas