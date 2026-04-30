"""
Módulo de Física do Simulador de Incêndio
Responsável por calcular o desenvolvimento do fogo
"""

# Tabela de constantes de crescimento do fogo (alfa) baseada na NFPA [kW/s^2]
ALFA_NFPA = {
    "lento": 0.00293,
    "medio": 0.01172,
    "rapido": 0.0469,
    "ultra_rapido": 0.1876
}

def calcular_hrr(tempo_s, crescimento="medio", limite_max_kw=None):
    """
    Calcula a Taxa de Liberação de Calor (Heat Release Rate - HRR) num dado instante t.
    
    Parâmetros:
    - tempo_s: O tempo atual da simulação em segundos.
    - crescimento: O tipo de fogo ('lento', 'medio', 'rapido', 'ultra_rapido').
    - limite_max_kw: Potência máxima que o fogo pode atingir (opcional).
    
    Retorna:
    - Q_dot: A potência do fogo em kW naquele segundo.
    """
    
    # Busca o valor de alfa no dicionário. Se o usuário digitar errado, usa 'medio' como padrão.
    alfa = ALFA_NFPA.get(crescimento, ALFA_NFPA["medio"])
    
    # A equação principal: Q_dot = alfa * t^2
    q_dot = alfa * (tempo_s ** 2)
    
    # O fogo não cresce ao infinito. Se a sala for pequena, ele estabiliza 
    # (limitação por falta de oxigênio ou porque o objeto já está todo em chamas).
    if limite_max_kw is not None and q_dot > limite_max_kw:
        return limite_max_kw
        
    return q_dot