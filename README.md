# App Web — Tubulão Térmico

Aplicação web para **curva analítica** e **regressão** a partir de dados de monitoramento em CSV, conforme a ideia do Gemini (arquitetura Frontend + Backend serverless/local).

## O que faz

- **Frontend:** upload de CSV (tempo em h, temperatura em °C), envio dos dados para a API e exibição do gráfico com dados experimentais e curva analítica ajustada.
- **Backend:** recebe tempos e temperaturas, executa a regressão (modelo Hill + solução analítica por Laplace/Stehfest) e devolve parâmetros otimizados e pontos para o gráfico.

## Estrutura

```
tubulaotermico_appweb/
├── README.md
├── requirements.txt          # servidor local (Flask, numpy, scipy)
├── server_local.py           # servidor Flask (frontend + /api/otimizar)
├── scripts/
│   └── test_cloud_function.py # teste da API (local ou Cloud Function)
├── static/
│   ├── index.html
│   └── app.js
└── backend/
    ├── main.py               # lógica de regressão + entrada Cloud Function (CORS)
    └── requirements.txt      # para deploy Google Cloud Functions
```

## Uso local

1. Instale as dependências:
   ```bash
   cd /home/flavio/Documentos/FEM/FLEXPDE/tubulaotermico_appweb
   pip install -r requirements.txt
   ```

2. Inicie o servidor:
   ```bash
   python server_local.py
   ```

3. Abra no navegador: **http://localhost:5000**

4. Selecione um CSV com colunas de tempo (h) e temperatura (°C), depois clique em **Calcular**. O gráfico e os parâmetros aparecem na mesma página.

## Formato do CSV

- Primeira coluna: tempo em horas (decimal com vírgula ou ponto).
- Segunda coluna: temperatura em °C.
- Cabeçalho opcional (ex.: `tempo (h),Temperatura (oC)`).
- Exemplo: linhas como `"0,01",25` ou `0.5,26.3` são aceitas.

## Deploy na nuvem (Google Cloud Functions)

O código em `backend/main.py` está preparado para **Google Cloud Functions** (Gen2), com **CORS** para o frontend em outro domínio.

### 1. Deploy via terminal (gcloud CLI)

Na pasta do projeto (`tubulaotermico_appweb/`):

```bash
gcloud functions deploy otimizar_tubulao \
  --gen2 \
  --runtime=python310 \
  --region=us-central1 \
  --source=./backend \
  --entry-point=otimizar_tubulao \
  --trigger-http \
  --allow-unauthenticated \
  --memory=512Mi \
  --timeout=60s
```

Após o deploy, anote a **URL da função** (ex.: `https://us-central1-PROJECT.cloudfunctions.net/otimizar_tubulao`).

### 2. Conectar o frontend à nuvem

Em `static/app.js`, defina a URL da Cloud Function:

```javascript
const API_OTIMIZAR_CLOUD = 'https://us-central1-SEU_PROJETO.cloudfunctions.net/otimizar_tubulao';
```

Deixe vazio para continuar usando o servidor local.

### 3. Testar a função na nuvem

```bash
# Testar URL local
python scripts/test_cloud_function.py

# Testar Cloud Function (substitua pela sua URL)
python scripts/test_cloud_function.py https://us-central1-PROJECT.cloudfunctions.net/otimizar_tubulao

# Ou via variável de ambiente
TEST_GCF_URL=https://... python scripts/test_cloud_function.py
```

### Dependências do backend (GCF)

O `backend/requirements.txt` já inclui: `functions-framework`, `numpy`, `scipy`. Não é necessário `pandas` para o cálculo.

### API

- **Corpo (POST):** `{ "tempos": [ ... ], "temperaturas": [ ... ], "chute": [ ... ] }` (chute opcional)
- **Resposta:** `parametros` (lista com nome, estimado, se, ic_inf, ic_sup, cv), `t_plot`, `T_plot`, `CI_lwr`, `CI_upr`, `erro_mae`

## Referência

Ideia e arquitetura descritas em:  
`tubulaotermico_output/TubulaoTermico/appweb_tubulaotermico_prompt.md`

## ALTERAÇÃO IMPORTANTE

Essa é, sem dúvida, a versão mais elegante e matematicamente madura do modelo.

Ao substituirmos as grandezas absolutas por **Coeficientes Adiabáticos ()** e a **Condutividade Relativa ()**, eliminamos completamente o ruído dimensional da equação e tornamos a regressão extremamente focada.

Para implementar exatamente o que conversamos:

1. **Reduzimos para 9 parâmetros**: `[dT_adi1, dT_adi2, tau1, beta1, tau2, beta2, k_rel, alpha1, alpha2]`.
2. **O consumo de cimento e o calor específico desaparecem** da formulação. Eles já estão embutidos nos $dT_{adi}$ (que agora são medidos diretamente em °C).
3. **Regressão em 2 Passos**: Usamos o pico do dado experimental para fatiar a curva. O Passo 1 ajusta só o modelo de Hill assumindo o tubulão como um meio infinito ($krel​=1$). O Passo 2 liga o "solo" e ajusta tudo.
4. **Derivada Exata:** Adicionei a função que calcula a velocidade de aquecimento ∂t∂T no domínio de Laplace (vˉ(s)=s⋅θˉ(s)) e devolve isso pro frontend.


