"""
Backend App Web Tubulão Térmico - MVP 9 Parâmetros (2 Passos + k_rel)
"""
import os
import math
import json
import numpy as np
import scipy.integrate as integrate
import scipy.stats as stats
from scipy.special import ive, kve
from scipy.optimize import least_squares

try:
    import functions_framework
except ImportError:
    functions_framework = None

DEFAULT_A = 0.45
DEFAULT_T_INI = 25.0

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

# Para Gauss-Laguerre (60 pontos) para integração rápida da Transformada de Laplace
_NODES_GL, _WEIGHTS_GL = np.polynomial.laguerre.laggauss(60)

def T_adi_hill(t, dT1, dT2, t1, b1, t2, b2):
    """Elevação Adiabática de Temperatura (°C) - Versão Vetorizada"""
    t_safe = np.where(t > 0, t, 1e-9)
    term1 = dT1 * (t_safe**b1) / (t1**b1 + t_safe**b1)
    term2 = dT2 * (t_safe**b2) / (t2**b2 + t_safe**b2)
    return np.where(t > 0, term1 + term2, 0.0)

def get_theta_bar_centro(s, params, a):
    """Retorna a transformada de Laplace do ganho de temperatura no centro.
    VETORIZADA sobre s."""
    dT1, dT2, t1, b1, t2, b2, k_rel, beta_alpha1, beta_alpha2 = params
    alpha1 = beta_alpha1 * 1e-4
    alpha2 = beta_alpha2 * 1e-4

    # Transformada de Laplace do calor adiabático via Gauss-Laguerre (Vetorizada sobre s)
    # integral = sum( w_j * f(x_j/s) / s )
    s_arr = np.atleast_1d(s)
    nodes_s = _NODES_GL.reshape(-1, 1) / s_arr.reshape(1, -1)
    T_nodes = T_adi_hill(nodes_s, dT1, dT2, t1, b1, t2, b2) 
    dT_adi_bar_s = np.sum(_WEIGHTS_GL.reshape(-1, 1) * T_nodes / s_arr.reshape(1, -1), axis=0)

    q1 = np.sqrt(s / alpha1)
    q2 = np.sqrt(s / alpha2)
    ratio_I = ive(1, q1 * a) / ive(0, q1 * a)
    ratio_K = kve(0, q2 * a) / kve(1, q2 * a)
    
    flux_ratio = k_rel * np.sqrt(alpha2 / alpha1)
    D_scaled = 1 + flux_ratio * ratio_I * ratio_K
    I0_a_scaled = ive(0, q1 * a)
    
    # Proteção contra divisões por zero ou NaNs em s muito grandes/pequenos
    mask_valid = (I0_a_scaled != 0) & (D_scaled != 0) & (~np.isinf(I0_a_scaled))
    term_sub = np.zeros_like(s)
    if isinstance(s, np.ndarray):
        term_sub[mask_valid] = (1.0 / I0_a_scaled[mask_valid]) * np.exp(-q1[mask_valid] * a) / D_scaled[mask_valid]
    else:
        if mask_valid:
            term_sub = (1.0 / I0_a_scaled) * np.exp(-q1 * a) / D_scaled
        
    return dT_adi_bar_s * (1 - term_sub)

def calc_temperatura_centro(tempos, params, T_ini=None, a=None):
    _T_ini = T_ini if T_ini is not None else DEFAULT_T_INI
    _a = a if a is not None else DEFAULT_A
    tempos = np.atleast_1d(tempos)
    res_T = np.full(tempos.shape, _T_ini, dtype=float)
    
    mask = tempos > 0
    t_val = tempos[mask]
    if len(t_val) == 0: return res_T

    ln2_t = np.log(2) / t_val
    i_indices = np.arange(1, 11).reshape(10, 1)
    S = i_indices * ln2_t.reshape(1, -1) # (10, N_t)
    
    theta_bar = get_theta_bar_centro(S, params, _a) # (10, N_t)
    soma = np.sum(V_STEHFEST.reshape(10, 1) * theta_bar, axis=0)
    res_T[mask] = _T_ini + soma * ln2_t
    return res_T

def calc_derivada_centro(tempos, params, a=None):
    _a = a if a is not None else DEFAULT_A
    tempos = np.atleast_1d(tempos)
    res_v = np.zeros(tempos.shape, dtype=float)
    
    mask = tempos > 0
    t_val = tempos[mask]
    if len(t_val) == 0: return res_v

    ln2_t = np.log(2) / t_val
    i_indices = np.arange(1, 11).reshape(10, 1)
    S = i_indices * ln2_t.reshape(1, -1)
    
    theta_bar = get_theta_bar_centro(S, params, _a)
    soma = np.sum(V_STEHFEST.reshape(10, 1) * S * theta_bar, axis=0)
    res_v[mask] = soma * ln2_t
    return res_v

# 9 Parâmetros: dT_adi1, dT_adi2, tau1, beta1, tau2, beta2, k_rel, beta_alpha1, beta_alpha2
# beta_alpha = alpha_real × 1e4  (ex: alpha=0.004 m²/h → beta=40)
DEFAULT_CHUTE = [45.0, 40.0, 10.0, 3.0, 25.0, 1.5, 2.9, 40.0, 30.0]
DEFAULT_BOUNDS_INF = [5.0,  5.0,  1.0, 0.5,  1.0, 0.5, 0.5,   5.0,  3.0]
DEFAULT_BOUNDS_SUP = [90.0, 90.0, 500.0, 10.0, 500.0, 10.0, 10.0, 100.0, 50.0]
N_PARAMS = 9
ALPHA_SCALE = 1e-4  # fator de conversão beta → alpha

# #region agent log
_DEBUG_LOG_PATH = os.environ.get("DEBUG_LOG_PATH", "debug.log")

def _to_jsonable(x):
    if hasattr(x, "tolist"):
        return [_to_jsonable(v) for v in x]
    if isinstance(x, (list, tuple)):
        return [_to_jsonable(v) for v in x]
    if isinstance(x, (np.floating, np.integer)):
        return float(x)
    return x

def _debug_log(msg, data):
    try:
        out = {"message": msg, "data": {k: _to_jsonable(v) for k, v in data.items()}}
        with open(_DEBUG_LOG_PATH, "a") as f:
            f.write(json.dumps(out) + "\n")
    except Exception:
        pass
# #endregion

def _to_beta_scale(vals):
    """Auto-converte alpha (índices 7,8) de unidades físicas para escala beta.
    Se o valor é pequeno (< 1.0), assume unidade física e multiplica por 1/ALPHA_SCALE."""
    vals = list(vals)
    for idx in (7, 8):
        if idx < len(vals) and vals[idx] < 1.0:
            vals[idx] /= ALPHA_SCALE  # ex: 0.004 → 40
    return vals

def run_otimizacao(tempos, temperaturas, chute=None, config=None):
    cfg = config or {}
    T_ini = float(cfg.get("T_ini", DEFAULT_T_INI))
    a = float(cfg["diametro"]) / 2.0 if "diametro" in cfg else float(cfg.get("raio", DEFAULT_A))
    chute = chute if chute is not None else list(cfg.get("chute", DEFAULT_CHUTE))
    bounds_inf = list(cfg.get("bounds_inf", DEFAULT_BOUNDS_INF))
    bounds_sup = list(cfg.get("bounds_sup", DEFAULT_BOUNDS_SUP))

    # Auto-converter de unidades físicas para escala beta se necessário
    chute = _to_beta_scale(chute)
    bounds_inf = _to_beta_scale(bounds_inf)
    bounds_sup = _to_beta_scale(bounds_sup)

    print(f"[DIAG] chute (beta):      {[f'{v:.2f}' for v in chute]}")
    print(f"[DIAG] bounds_inf (beta):  {[f'{v:.2f}' for v in bounds_inf]}")
    print(f"[DIAG] bounds_sup (beta):  {[f'{v:.2f}' for v in bounds_sup]}")
    print(f"[DIAG] T_ini={T_ini}, a={a}")

    t_exp = np.asarray(tempos, dtype=float)
    T_exp = np.asarray(temperaturas, dtype=float)
    valid = (t_exp > 0.1) & (~np.isnan(t_exp)) & (~np.isnan(T_exp))
    t_exp, T_exp = t_exp[valid], T_exp[valid]

    idx_pico_exp = np.argmax(T_exp)
    idx_fit = np.unique(np.round(np.concatenate([
        np.linspace(0, idx_pico_exp, 40), np.linspace(idx_pico_exp + 1, len(t_exp) - 1, 60)
    ])).astype(int))
    idx_fit = idx_fit[(idx_fit >= 0) & (idx_fit < len(t_exp))]
    t_fit, T_fit = t_exp[idx_fit], T_exp[idx_fit]

    # ==========================================================
    # PASSO 1: Regressão na Fase de Aquecimento (Meio Infinito)
    # ==========================================================
    idx_step1 = [i for i in idx_fit if i <= idx_pico_exp + 5] # Até o pico + pequena margem
    t_step1, T_step1 = t_exp[idx_step1], T_exp[idx_step1]

    def residuals_step1(p_hill):
        # Fixa k_rel = 1.0 e alpha1 = alpha2 (fronteira térmica invisível)
        p_full = list(p_hill) + [1.0, chute[7], chute[7]]
        return calc_temperatura_centro(t_step1, p_full, T_ini=T_ini, a=a) - T_step1

    chute_step1 = np.clip(chute[:6], bounds_inf[:6], bounds_sup[:6])
    print("\n[OTIMIZAÇÃO] Iniciando Passo 1 (Aquecimento)...")
    res1 = least_squares(residuals_step1, chute_step1, bounds=(bounds_inf[:6], bounds_sup[:6]),
                         method="trf", x_scale="jac", diff_step=1e-6, ftol=1e-5, xtol=1e-5, verbose=2)
    p_hill_opt = res1.x

    # ==========================================================
    # PASSO 2: Regressão Completa (Solo ativado)
    # ==========================================================
    chute_step2 = np.array(list(p_hill_opt) + [chute[6], chute[7], chute[8]])
    chute_step2 = np.clip(chute_step2, bounds_inf, bounds_sup)
    _debug_log("Passo2 entrada", {"len_chute_step2": len(chute_step2), "len_bounds_inf": len(bounds_inf), "len_bounds_sup": len(bounds_sup), "chute_step2": chute_step2})
    def residuals_step2(p_full):
        return calc_temperatura_centro(t_fit, p_full, T_ini=T_ini, a=a) - T_fit

    print("\n[OTIMIZAÇÃO] Iniciando Passo 2 (Completo 9D)...")
    res2 = least_squares(residuals_step2, chute_step2, bounds=(bounds_inf, bounds_sup),
                         method="trf", x_scale="jac", diff_step=1e-6, ftol=1e-5, xtol=1e-5, verbose=2)
    p_opt = res2.x
    # Verificação cruzada: custo real vs custo reportado
    resid_check = calc_temperatura_centro(t_fit, p_opt, T_ini=T_ini, a=a) - T_fit
    cost_check = 0.5 * np.sum(resid_check**2)
    igual_678 = [float(p_opt[j]) == float(chute_step2[j]) for j in (6, 7, 8)]
    print(f"[DIAG] Step1 result:  {[f'{v:.2f}' for v in p_hill_opt]}")
    print(f"[DIAG] Step2 chute:   {[f'{v:.2f}' for v in chute_step2]}")
    print(f"[DIAG] Step2 result:  {[f'{v:.2f}' for v in p_opt]}")
    print(f"[DIAG] Step2 nfev={res2.nfev}, cost_scipy={res2.cost:.4e}, cost_check={cost_check:.4e}, status={res2.status}")
    print(f"[DIAG] k_rel/a1/a2 unchanged: {igual_678}")
    _debug_log("Passo2 saida", {"len_p_opt": len(p_opt), "p_opt": p_opt, "chute_step2": chute_step2,
                                  "p_opt_igual_chute_em_678": igual_678, "res2_nfev": res2.nfev,
                                  "res2_cost": res2.cost, "res2_optimality": res2.optimality,
                                  "res2_status": res2.status, "res2_message": res2.message})

    # --- Estatística Assintótica usando Jacobiano do scipy ---
    n_obs, p_par = len(t_fit), len(p_opt)
    df_resid = n_obs - p_par
    t_crit = stats.t.ppf(1 - 0.05 / 2, df_resid)

    # Jacobiano do scipy + residuais verificados (cross-check)
    Fdot = res2.jac

    _debug_log("Jacobiano stats", {
        "Fdot_col_norms": np.linalg.norm(Fdot, axis=0).tolist(),
        "Fdot_shape": list(Fdot.shape),
    })

    # Usar resíduos verificados (avaliados diretamente), não res2.fun
    # que pode estar em escala interna diferente se x_scale foi usado
    s2 = np.sum(resid_check**2) / df_resid
    s = np.sqrt(s2)

    scale_factors = np.linalg.norm(Fdot, axis=0)
    scale_factors[scale_factors < 1e-12] = 1.0
    F_scaled = Fdot / scale_factors
    FtF_s = F_scaled.T @ F_scaled

    FtF_s_inv = np.linalg.pinv(FtF_s, rcond=1e-5)
    D_inv = np.diag(1.0 / scale_factors)
    FtF_inv = D_inv @ FtF_s_inv @ D_inv

    Sigma = s2 * FtF_inv
    SE_param = np.sqrt(np.clip(np.diag(Sigma), 0, None))
    IC_inf, IC_sup = p_opt - t_crit * SE_param, p_opt + t_crit * SE_param

    CV_pct = np.zeros(p_par)
    for idx_cv in range(p_par):
        if abs(p_opt[idx_cv]) > 1e-10:
            CV_pct[idx_cv] = (SE_param[idx_cv] / abs(p_opt[idx_cv])) * 100

    # --- Amostragem Final ---
    indices_plot = np.unique(np.concatenate([np.where(t_exp < 2.0)[0], np.linspace(0, idx_pico_exp, 50, dtype=int), np.linspace(idx_pico_exp, len(t_exp)-1, 80, dtype=int)]))
    indices_plot = indices_plot[indices_plot < len(t_exp)]
    t_plot = t_exp[indices_plot]
    T_plot = calc_temperatura_centro(t_plot, p_opt, T_ini=T_ini, a=a)
    v_plot = calc_derivada_centro(t_plot, p_opt, a=a)

    # Bandas usando eps_rel maior para diferenças fintas da curva (separado da estatística)
    eps_band = 1e-2
    Fdot_grade = np.zeros((len(t_plot), p_par))
    for j in range(p_par):
        h = max(eps_band * abs(p_opt[j]), 1e-4)
        p_plus, p_minus = p_opt.copy(), p_opt.copy()
        p_plus[j] += h
        p_minus[j] -= h
        T_plus = calc_temperatura_centro(t_plot, p_plus, T_ini=T_ini, a=a)
        T_minus = calc_temperatura_centro(t_plot, p_minus, T_ini=T_ini, a=a)
        Fdot_grade[:, j] = (T_plus - T_minus) / (2 * h)

    se_curva = s * np.sqrt(np.clip(np.sum((Fdot_grade @ FtF_inv) * Fdot_grade, axis=1), 0, None))
    CI_lwr, CI_upr = T_plot - t_crit * se_curva, T_plot + t_crit * se_curva

    nomes_parametros = ["dT_adi1", "dT_adi2", "tau1", "beta1", "tau2", "beta2", "k_rel", "alpha1", "alpha2"]
    stats_data = []
    for i in range(p_par):
        est, ic_i, ic_s, se = float(p_opt[i]), float(IC_inf[i]), float(IC_sup[i]), float(SE_param[i])
        # Converter beta → alpha (unidades físicas) para exibição
        if i >= 7:
            est *= ALPHA_SCALE
            ic_i *= ALPHA_SCALE
            ic_s *= ALPHA_SCALE
            se *= ALPHA_SCALE
        stats_data.append({"nome": nomes_parametros[i], "estimado": est, "se": se, "ic_inf": ic_i, "ic_sup": ic_s, "cv": float(CV_pct[i])})

    return {
        "parametros": stats_data,
        "t_plot": t_plot.tolist(),
        "T_plot": T_plot.tolist(),
        "v_plot": v_plot.tolist(),
        "CI_lwr": CI_lwr.tolist(),
        "CI_upr": CI_upr.tolist(),
        "erro_mae": float(np.mean(np.abs(res2.fun)))
    }


def run_curva(params, config=None, tempos=None):
    """Gera apenas a curva T x t (9 parâmetros). Sem regressão."""
    params = np.asarray(params, dtype=float).copy()
    if len(params) != N_PARAMS:
        return {"error": "params deve ter 9 elementos."}
    # Auto-converter alpha de unidades físicas para escala beta
    for idx in (7, 8):
        if params[idx] < 1.0:
            params[idx] /= ALPHA_SCALE

    cfg = config or {}
    T_ini = float(cfg.get("T_ini", DEFAULT_T_INI))
    a = float(cfg["diametro"]) / 2.0 if "diametro" in cfg else float(cfg.get("raio", DEFAULT_A))
    if tempos is None or len(tempos) == 0:
        tempos = np.linspace(0.1, 100.0, 300)
    else:
        tempos = np.asarray(tempos, dtype=float)
    T_plot = calc_temperatura_centro(tempos, params, T_ini=T_ini, a=a)
    return {"t_plot": tempos.tolist(), "T_plot": T_plot.tolist()}


if functions_framework is not None:
    @functions_framework.http
    def otimizar_tubulao(request):
        if request.method == "OPTIONS":
            return ("", 204, {"Access-Control-Allow-Origin": "*", "Access-Control-Allow-Methods": "POST", "Access-Control-Allow-Headers": "Content-Type"})
        req = request.get_json(silent=True)
        return (run_otimizacao(req.get("tempos", []), req.get("temperaturas", []), req.get("chute"), req.get("config")), 200, {"Access-Control-Allow-Origin": "*"})