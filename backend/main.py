"""
Backend do App Web Tubulão Térmico - MVP 10 Parâmetros + Estatística.
"""
import math
import numpy as np
import scipy.integrate as integrate
import scipy.stats as stats
from scipy.special import ive, kve
from scipy.optimize import least_squares

try:
    import functions_framework
except ImportError:
    functions_framework = None

# --- Constantes padrão do modelo (podem ser sobrescritas via config no request) ---
DEFAULT_A = 0.45
DEFAULT_T_INI = 25.0
DEFAULT_C_CIM = 300

def get_stehfest_V(n=10):
    V = np.zeros(n)
    for i in range(1, n + 1):
        k_min = (i + 1) // 2
        k_max = min(i, n // 2)
        sum_k = 0
        for k in range(k_min, k_max + 1):
            term = ((k ** (n // 2) * math.factorial(2 * k)) /
                    (math.factorial(n // 2 - k) * math.factorial(k) * math.factorial(k - 1) * math.factorial(i - k) * math.factorial(2 * k - i)))
            sum_k += term
        V[i - 1] = ((-1) ** (i + n // 2)) * sum_k
    return V

V_STEHFEST = get_stehfest_V(10)

def Q_t_hill(t, Qi1, Qi2, t1, b1, t2, b2, C_cim=None):
    C = C_cim if C_cim is not None else DEFAULT_C_CIM
    if t <= 0:
        return 0.0
    term1 = Qi1 * (t**b1) / (t1**b1 + t**b1)
    term2 = Qi2 * (t**b2) / (t2**b2 + t**b2)
    return C * (term1 + term2)


def calc_temperatura_centro(tempos, params, T_ini=None, a=None, C_cim=None):
    """Curva analítica. params com 10 elementos. T_ini, a, C_cim opcionais (fixos)."""
    _T_ini = T_ini if T_ini is not None else DEFAULT_T_INI
    _a = a if a is not None else DEFAULT_A
    _C_cim = C_cim if C_cim is not None else DEFAULT_C_CIM

    Qi1, Qi2, t1, b1, t2, b2, k1, alpha1, k2, alpha2 = params
    rho_cp1 = k1 / alpha1

    def get_g_bar(s):
        lim_sup = min(2000 / s, 1e5)
        res, _ = integrate.quad(
            lambda t: Q_t_hill(t, Qi1, Qi2, t1, b1, t2, b2, C_cim=_C_cim) * np.exp(-s * t),
            0, lim_sup, limit=100,
        )
        return s * res

    def theta_bar_centro(s):
        q1 = np.sqrt(s / alpha1)
        q2 = np.sqrt(s / alpha2)
        g_bar_s = get_g_bar(s)
        ratio_I = ive(1, q1 * _a) / ive(0, q1 * _a)
        ratio_K = kve(0, q2 * _a) / kve(1, q2 * _a)
        D_scaled = 1 + (k1 * q1) / (k2 * q2) * ratio_I * ratio_K
        I0_a_scaled = ive(0, q1 * _a)
        if np.isinf(I0_a_scaled) or np.isnan(I0_a_scaled) or D_scaled == 0:
            term_sub = 0
        else:
            term_sub = (1.0 / I0_a_scaled) * np.exp(-q1 * _a) / D_scaled
        return (g_bar_s / (s * rho_cp1)) * (1 - term_sub)

    res_T = np.zeros(len(tempos))
    for j, t in enumerate(tempos):
        if t <= 0:
            res_T[j] = _T_ini
            continue
        ln2_t = np.log(2) / t
        soma = sum(V_STEHFEST[i - 1] * theta_bar_centro(i * ln2_t) for i in range(1, 11))
        res_T[j] = _T_ini + soma * ln2_t
    return res_T

# Valores padrão para chute e limites (10 parâmetros; C_cim é fixo via config)
DEFAULT_CHUTE = [113.16, 93.84, 15.0, 1.9, 250.0, 1.5, 9.36, 0.0045, 3.0, 0.002]
DEFAULT_BOUNDS_INF = [50.0, 50.0, 1.0, 0.5, 10.0, 0.5, 2.0, 0.001, 0.5, 0.0005]
DEFAULT_BOUNDS_SUP = [200.0, 200.0, 50.0, 5.0, 500.0, 5.0, 15.0, 0.015, 8.0, 0.010]
N_PARAMS = 10


def run_otimizacao(tempos, temperaturas, chute=None, config=None):
    cfg = config or {}
    T_ini = float(cfg.get("T_ini", DEFAULT_T_INI))
    if "diametro" in cfg:
        a = float(cfg["diametro"]) / 2.0
    else:
        a = float(cfg.get("raio", DEFAULT_A))
    C_cim = float(cfg.get("C_cim", DEFAULT_C_CIM))
    chute = chute if chute is not None else list(cfg.get("chute", DEFAULT_CHUTE))
    bounds_inf = list(cfg.get("bounds_inf", DEFAULT_BOUNDS_INF))
    bounds_sup = list(cfg.get("bounds_sup", DEFAULT_BOUNDS_SUP))
    confianca = float(cfg.get("confianca", 0.95))
    eps_rel = float(cfg.get("eps_rel", 1e-4))

    t_exp = np.asarray(tempos, dtype=float)
    T_exp = np.asarray(temperaturas, dtype=float)

    valid = (t_exp > 0.1) & (~np.isnan(t_exp)) & (~np.isnan(T_exp))
    t_exp, T_exp = t_exp[valid], T_exp[valid]

    if len(t_exp) < 15:
        return {"error": "Dados insuficientes (mínimo 15 pontos válidos)."}

    if len(chute) != N_PARAMS or len(bounds_inf) != N_PARAMS or len(bounds_sup) != N_PARAMS:
        return {"error": "chute e bounds devem ter 10 elementos cada."}

    chute = np.asarray(chute, dtype=float)
    bounds_inf = np.asarray(bounds_inf, dtype=float)
    bounds_sup = np.asarray(bounds_sup, dtype=float)

    idx_pico = np.argmax(T_exp)
    idx_fit = np.unique(np.round(np.concatenate([
        np.linspace(0, idx_pico, 40), np.linspace(idx_pico + 1, len(t_exp) - 1, 60)
    ])).astype(int))
    idx_fit = idx_fit[(idx_fit >= 0) & (idx_fit < len(t_exp))]
    t_fit, T_fit = t_exp[idx_fit], T_exp[idx_fit]

    def residuals(p):
        return calc_temperatura_centro(t_fit, p, T_ini=T_ini, a=a, C_cim=C_cim) - T_fit

    res = least_squares(residuals, chute, bounds=(bounds_inf, bounds_sup), method="trf")
    p_opt = res.x

    # --- Estatística Assintótica ---
    n_obs = len(t_fit)
    p_par = len(p_opt)
    df_resid = n_obs - p_par
    alpha = 1.0 - confianca
    t_crit = stats.t.ppf(1 - alpha / 2, df_resid)

    Fdot = np.zeros((n_obs, p_par))
    for j in range(p_par):
        h = max(eps_rel * abs(p_opt[j]), 1e-8)
        p_plus, p_minus = p_opt.copy(), p_opt.copy()
        p_plus[j] += h; p_minus[j] -= h
        Fdot[:, j] = (
            calc_temperatura_centro(t_fit, p_plus, T_ini=T_ini, a=a, C_cim=C_cim)
            - calc_temperatura_centro(t_fit, p_minus, T_ini=T_ini, a=a, C_cim=C_cim)
        ) / (2 * h)

    residuos_opt = res.fun
    SSR = np.sum(residuos_opt**2)
    s2 = SSR / df_resid
    s = np.sqrt(s2)

    # =================================================================
    # CORREÇÃO: Pré-condicionamento (Scaling) do Jacobiano
    # =================================================================
    # 1. Calcula a norma (tamanho) de cada coluna do Jacobiano
    scale_factors = np.linalg.norm(Fdot, axis=0)
    # Evita divisão por zero se alguma derivada for nula
    scale_factors[scale_factors < 1e-12] = 1.0

    # 2. Divide F pelo fator de escala (todas as colunas passam a ter norma 1)
    F_scaled = Fdot / scale_factors

    # 3. Matriz de correlação (bem condicionada)
    FtF_s = F_scaled.T @ F_scaled

    try:
        cond_FtF = np.linalg.cond(FtF_s)
        if not np.isfinite(cond_FtF) or cond_FtF > 1e30:
            cond_FtF = 1e30
    except Exception:
        cond_FtF = 1e30

    # 4. Inversão robusta usando SVD (Pseudo-inversa truncada)
    # O rcond=1e-5 joga fora direções de multicolinearidade perfeita que fazem o IC explodir
    FtF_s_inv = np.linalg.pinv(FtF_s, rcond=1e-5)

    # 5. Desfaz o escalonamento para voltar à dimensão real dos parâmetros
    # (F'F)^-1 = D^-1 * (F_s'F_s)^-1 * D^-1
    D_inv = np.diag(1.0 / scale_factors)
    FtF_inv = D_inv @ FtF_s_inv @ D_inv
    # =================================================================

    Sigma = s2 * FtF_inv
    diag_Sigma = np.diag(Sigma)

    # Garante que variâncias não sejam negativas por erro de arredondamento
    diag_Sigma = np.clip(diag_Sigma, 0, None)
    SE_param = np.sqrt(diag_Sigma)

    IC_inf = p_opt - t_crit * SE_param
    IC_sup = p_opt + t_crit * SE_param

    # Cálculo seguro do Coeficiente de Variação
    CV_pct = np.zeros(p_par)
    for idx_cv in range(p_par):
        if abs(p_opt[idx_cv]) > 1e-10:
            CV_pct[idx_cv] = (SE_param[idx_cv] / abs(p_opt[idx_cv])) * 100
        else:
            CV_pct[idx_cv] = 0.0

    CV_pct = np.clip(CV_pct, 0, 1e6)

    # --- Amostragem e Bandas da Curva ---
    indices_plot = np.unique(np.concatenate([
        np.where(t_exp < 2.0)[0],
        np.linspace(0, idx_pico, 50, dtype=int),
        np.linspace(idx_pico, len(t_exp) - 1, 80, dtype=int)
    ]))
    indices_plot = indices_plot[indices_plot < len(t_exp)]
    t_plot = t_exp[indices_plot]
    T_plot = calc_temperatura_centro(t_plot, p_opt, T_ini=T_ini, a=a, C_cim=C_cim)

    Fdot_grade = np.zeros((len(t_plot), p_par))
    for j in range(p_par):
        h = max(eps_rel * abs(p_opt[j]), 1e-8)
        p_plus, p_minus = p_opt.copy(), p_opt.copy()
        p_plus[j] += h; p_minus[j] -= h
        Fdot_grade[:, j] = (
            calc_temperatura_centro(t_plot, p_plus, T_ini=T_ini, a=a, C_cim=C_cim)
            - calc_temperatura_centro(t_plot, p_minus, T_ini=T_ini, a=a, C_cim=C_cim)
        ) / (2 * h)

    var_curva = np.sum((Fdot_grade @ FtF_inv) * Fdot_grade, axis=1)
    var_curva = np.clip(var_curva, 0, None)
    se_curva = s * np.sqrt(var_curva)
    se_curva = np.where(np.isfinite(se_curva), se_curva, 0.0)
    se_curva = np.clip(se_curva, 0, 1e10)
    CI_lwr = T_plot - t_crit * se_curva
    CI_upr = T_plot + t_crit * se_curva

    nomes_parametros = ["Qi1", "Qi2", "tau1", "beta1", "tau2", "beta2", "k1", "alpha1", "k2", "alpha2"]

    stats_aviso = None
    if cond_FtF > 1e10:
        stats_aviso = f"Matriz F'F mal condicionada (cond={cond_FtF:.2e}). IC e bandas podem ser pouco confiáveis."
    elif cond_FtF > 1e6:
        stats_aviso = f"Matriz F'F com condicionamento alto (cond={cond_FtF:.2e}). Use os IC com cautela."

    # Preparando dados estatísticos para o frontend
    stats_data = []
    for i in range(p_par):
        stats_data.append({
            "nome": nomes_parametros[i],
            "estimado": float(p_opt[i]),
            "se": float(SE_param[i]),
            "ic_inf": float(IC_inf[i]),
            "ic_sup": float(IC_sup[i]),
            "cv": float(CV_pct[i])
        })

    return {
        "parametros": stats_data,
        "t_plot": t_plot.tolist(),
        "T_plot": T_plot.tolist(),
        "CI_lwr": CI_lwr.tolist(),
        "CI_upr": CI_upr.tolist(),
        "erro_mae": float(np.mean(np.abs(res.fun))),
        "confianca": confianca,
        "cond_FtF": float(cond_FtF),
        "stats_aviso": stats_aviso,
    }


def run_curva(params, config=None, tempos=None):
    """
    Gera apenas a curva temperatura x tempo (sem regressão).
    params: lista de 10 parâmetros [Qi1, Qi2, tau1, beta1, tau2, beta2, k1, alpha1, k2, alpha2].
    config: opcional { T_ini, diametro ou raio, C_cim }.
    tempos: opcional array de tempos (h). Se None, usa grade padrão (0.1 a 100 h, 300 pontos).
    """
    params = np.asarray(params, dtype=float)
    if len(params) != N_PARAMS:
        return {"error": "params deve ter 10 elementos."}

    cfg = config or {}
    T_ini = float(cfg.get("T_ini", DEFAULT_T_INI))
    if "diametro" in cfg:
        a = float(cfg["diametro"]) / 2.0
    else:
        a = float(cfg.get("raio", DEFAULT_A))
    C_cim = float(cfg.get("C_cim", DEFAULT_C_CIM))

    if tempos is None or len(tempos) == 0:
        tempos = np.linspace(0.1, 100.0, 300)
    else:
        tempos = np.asarray(tempos, dtype=float)

    T_plot = calc_temperatura_centro(tempos, params, T_ini=T_ini, a=a, C_cim=C_cim)
    return {
        "t_plot": tempos.tolist(),
        "T_plot": T_plot.tolist(),
    }


if functions_framework is not None:
    @functions_framework.http
    def otimizar_tubulao(request):
        # CORS: permite que o app web (localhost ou outro domínio) chame a API
        if request.method == "OPTIONS":
            headers = {
                "Access-Control-Allow-Origin": "*",
                "Access-Control-Allow-Methods": "POST, OPTIONS",
                "Access-Control-Allow-Headers": "Content-Type",
                "Access-Control-Max-Age": "3600",
            }
            return ("", 204, headers)

        headers = {"Access-Control-Allow-Origin": "*"}

        request_json = request.get_json(silent=True)
        if not request_json or "tempos" not in request_json:
            return ({"error": "Dados insuficientes."}, 400, headers)

        out = run_otimizacao(
            request_json.get("tempos", []),
            request_json.get("temperaturas", []),
            chute=request_json.get("chute"),
            config=request_json.get("config"),
        )
        if "error" in out:
            return (out, 400, headers)
        return (out, 200, headers)

    @functions_framework.http
    def curva_tubulao(request):
        """Gera apenas a curva T x t (sem regressão). Corpo: { params, config?, tempos? }."""
        if request.method == "OPTIONS":
            headers = {
                "Access-Control-Allow-Origin": "*",
                "Access-Control-Allow-Methods": "POST, OPTIONS",
                "Access-Control-Allow-Headers": "Content-Type",
                "Access-Control-Max-Age": "3600",
            }
            return ("", 204, headers)
        headers = {"Access-Control-Allow-Origin": "*"}
        request_json = request.get_json(silent=True)
        if not request_json or "params" not in request_json:
            return ({"error": "Envie 'params' (lista de 10 valores)."}, 400, headers)
        out = run_curva(
            request_json["params"],
            config=request_json.get("config"),
            tempos=request_json.get("tempos"),
        )
        if "error" in out:
            return (out, 400, headers)
        return (out, 200, headers)