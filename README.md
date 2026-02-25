# Tubulão Térmico - Análise de Aquecimento de Concreto

![Status: Stable](https://img.shields.io/badge/Status-Estável-green)
![Math: Laplace Transform](https://img.shields.io/badge/Matemática-Transformada_de_Laplace-blue)
![Backend: Python](https://img.shields.io/badge/Backend-Python-yellow)

Este projeto é uma ferramenta avançada para a simulação e análise do comportamento térmico de tubulões de concreto. Ele permite gerar curvas teóricas de aquecimento adiabático e realizar regressões matemáticas complexas a partir de dados experimentais de monitoramento.

## 🚀 Funcionalidades Principais

### 1. Monitoramento (Gerar Curva)
Simula a evolução da temperatura no centro do tubulão com base em parâmetros térmicos conhecidos.
- Visualização instantânea da curva de aquecimento.
- Cálculo preciso da velocidade de aquecimento ($\partial T / \partial t$).
- Suporte a modelos bi-fásicos (duas funções de Hill combinadas).

### 2. Análise (Fazer Regressão)
Ajusta automaticamente os parâmetros físicos do modelo a partir de um arquivo CSV de monitoramento real.
- Algoritmo de otimização em dois passos para máxima convergência.
- Cálculo de bandas de confiança e erros estatísticos.
- Exportação visual dos resultados comparando o dado real vs. modelo teórico.

## 🧠 Fundamentação Matemática

O "motor" deste app utiliza uma abordagem híbrida analítico-numérica para garantir velocidade e precisão:

- **Modelo Cinético**: Baseado em duas funções de Hill acopladas para descrever a liberação de calor do cimento.
- **Solução no Domínio de Laplace**: As equações diferenciais de condução de calor são resolvidas analiticamente no domínio de $s$ para considerar a geometria do tubulão e a interface com o solo.
- **Inversão de Stehfest**: Algoritmo robusto para retornar do domínio de Laplace para o domínio do tempo.
- **Integração de Gauss-Legendre**: Implementação otimizada (40x mais rápida que métodos tradicionais) para resolver a integral de convolução adiabática via quadratura numérica.

## 🛠️ Tecnologias Utilizadas

- **Backend**: Python 3.10 com Flask.
- **Processamento**: NumPy e SciPy (Otimização e álgebra linear).
- **Frontend**: JavaScript (Vanilla ES6).
- **Gráficos**: Chart.js com suporte a Zoom e Bandas de Confiança.

## 📂 Formato de Dados (CSV)

O sistema aceita arquivos CSV com as seguintes características:
- **Coluna 1**: Tempo em horas (aceita vírgula ou ponto, ex: `0,5` ou `0.5`).
- **Coluna 2**: Temperatura medida em °C.
- **Cabeçalho**: Opcional.

Exemplo de conteúdo:
```csv
tempo,temperatura
0,25.4
0.5,26.1
1.0,27.5
```

## 💻 Execução Local

Para rodar o projeto em sua máquina:

1. **Instale as dependências**:
   ```bash
   pip install -r requirements.txt
   ```

2. **Inicie o servidor**:
   ```bash
   python app.py
   ```

3. **Acesse no navegador**:
   `http://localhost:5000`

---
**Desenvolvido para análise avançada de estruturas de fundação em concreto de grande volume.**
