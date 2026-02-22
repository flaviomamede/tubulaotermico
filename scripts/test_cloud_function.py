#!/usr/bin/env python3
"""
Testa se a Cloud Function (ou o servidor local) está respondendo corretamente.
Uso:
  python scripts/test_cloud_function.py
  python scripts/test_cloud_function.py https://us-central1-PROJECT.cloudfunctions.net/otimizar_tubulao
  TEST_GCF_URL=https://... python scripts/test_cloud_function.py
"""
import json
import os
import sys
import urllib.request

# Dados mínimos para teste (15+ pontos, tempo > 0.1)
tempos = [0.2, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 9.0, 10.0, 12.0, 15.0, 20.0]
temperaturas = [25.1, 25.3, 25.8, 26.2, 26.5, 26.8, 27.0, 27.2, 27.3, 27.4, 27.5, 27.6, 27.7, 27.7, 27.8, 27.8, 27.9, 27.9, 28.0]

def main():
    url = (
        os.environ.get("TEST_GCF_URL")
        or (sys.argv[1] if len(sys.argv) > 1 else None)
        or "http://127.0.0.1:5000/api/otimizar"
    )
    payload = json.dumps({"tempos": tempos, "temperaturas": temperaturas}).encode("utf-8")
    req = urllib.request.Request(
        url,
        data=payload,
        method="POST",
        headers={"Content-Type": "application/json"},
    )
    try:
        with urllib.request.urlopen(req, timeout=120) as resp:
            body = resp.read().decode("utf-8")
            data = json.loads(body)
            if "error" in data:
                print("ERRO da API:", data["error"])
                sys.exit(1)
            print("OK Status:", resp.status)
            print("Parametros:", len(data.get("parametros", [])))
            print("erro_mae:", data.get("erro_mae"))
            print("t_plot length:", len(data.get("t_plot", [])))
            print("CI_lwr/CI_upr:", "OK" if data.get("CI_lwr") and data.get("CI_upr") else "ausente")
    except urllib.error.HTTPError as e:
        body = e.read().decode("utf-8")
        print("HTTP Error", e.code, body[:500])
        sys.exit(1)
    except Exception as e:
        print("Falha:", e)
        sys.exit(1)


if __name__ == "__main__":
    main()
