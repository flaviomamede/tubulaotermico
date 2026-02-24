"""
Servidor local para o App Web Tubulão Térmico.
Serve o frontend estático e a API de otimização (mesma lógica do backend/main.py).
Uso: python server_local.py  →  http://localhost:5000
"""
import json
import os
from flask import Flask, request, send_from_directory

from backend.main import run_otimizacao, run_curva

app = Flask(__name__, static_folder="static", static_url_path="")
app.config["SEND_FILE_MAX_AGE_DEFAULT"] = 0


@app.route("/")
def index():
    return send_from_directory(app.static_folder, "index.html")


@app.route("/<path:path>")
def static_files(path):
    return send_from_directory(app.static_folder, path)


@app.route("/api/otimizar", methods=["POST"])
def api_otimizar():
    """Recebe JSON com tempos e temperaturas; devolve parâmetros e pontos do gráfico."""
    if request.method != "POST":
        return {"error": "Método não permitido"}, 405
    try:
        data = request.get_json(force=True, silent=True) or {}
    except Exception:
        return {"error": "JSON inválido"}, 400

    tempos = data.get("tempos", [])
    temperaturas = data.get("temperaturas", [])
    chute = data.get("chute")
    config = data.get("config")

    if not tempos or not temperaturas:
        return {"error": "Envie 'tempos' e 'temperaturas' no corpo da requisição."}, 400

    # DUMP THE EXACT PAYLOAD FOR DIAGNOSIS
    with open("/home/flavio/Documentos/FEM/FLEXPDE/.cursor/payload_dump.json", "w") as f:
        json.dump(data, f)
    print(f"[HTTP] Dumped request with {len(tempos)} pts to payload_dump.json")

    out = run_otimizacao(tempos, temperaturas, chute=chute, config=config)
    if "error" in out:
        return out, 400
    return out


@app.route("/api/curva", methods=["POST"])
def api_curva():
    """Gera apenas a curva T x t. Corpo: { params, config?, tempos? }."""
    if request.method != "POST":
        return {"error": "Método não permitido"}, 405
    try:
        data = request.get_json(force=True, silent=True) or {}
    except Exception:
        return {"error": "JSON inválido"}, 400
    params = data.get("params", [])
    if not params:
        return {"error": "Envie 'params' (lista de 9 parâmetros)."}, 400
    out = run_curva(params, config=data.get("config"), tempos=data.get("tempos"))
    if "error" in out:
        return out, 400
    return out


if __name__ == "__main__":
    port = int(os.environ.get("PORT", 5000))
    app.run(host="0.0.0.0", port=port, debug=True)
