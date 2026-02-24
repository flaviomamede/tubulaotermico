import time
import math
import numpy as np
from scipy import integrate
from scipy.special import ive, kve
import matplotlib.pyplot as plt

# --- CONFIGURAÇÃO E CONSTANTES ---
V_STEHFEST = np.array([
    -0.0013888888888888889, 0.22361111111111112, -4.95, 34.0, -102.375, 
    163.8, -143.1, 68.4, -16.875, 1.6875
])

# Para Gauss-Legendre (150 pontos) no intervalo [0, 1]
NODES_GLeg, WEIGHTS_GLeg = np.polynomial.legendre.leggauss(150)
NODES_U = 0.5 * (NODES_GLeg + 1)
WEIGHTS_U = 0.5 * WEIGHTS_GLeg

def T_adi_hill_stable(t, dT1, dT2, t1, b1, t2, b2):
    """Fórmula estável: 1 / (1 + (tau/t)^b)"""
    with np.errstate(divide='ignore', invalid='ignore'):
        t_m = np.maximum(t, 1e-12)
        res1 = dT1 / (1.0 + (t1 / t_m)**b1)
        res2 = dT2 / (1.0 + (t2 / t_m)**b2)
        res = res1 + res2
    return np.where(t > 0, res, 0.0)

def T_adi_hill_OLD(t, dT1, dT2, t1, b1, t2, b2):
    if t <= 0: return 0.0
    return dT1 * (t**b1) / (t1**b1 + t**b1) + dT2 * (t**b2) / (t2**b2 + t**b2)

# --- VERSÃO ANTIGA (COM QUAD) ---
def get_theta_bar_centro_OLD(s, params, a):
    dT1, dT2, t1, b1, t2, b2, k_rel, beta_alpha1, beta_alpha2 = params
    alpha1 = beta_alpha1 * 1e-4
    alpha2 = beta_alpha2 * 1e-4
    lim_sup = min(2000 / s, 1e5)
    dT_adi_bar_s, _ = integrate.quad(
        lambda t: T_adi_hill_OLD(t, dT1, dT2, t1, b1, t2, b2) * np.exp(-s * t),
        0, lim_sup, limit=100
    )
    q1 = np.sqrt(s / alpha1)
    q2 = np.sqrt(s / alpha2)
    ratio_I = ive(1, q1 * a) / ive(0, q1 * a)
    ratio_K = kve(0, q2 * a) / kve(1, q2 * a)
    flux_ratio = k_rel * np.sqrt(alpha2 / alpha1)
    D_scaled = 1 + flux_ratio * ratio_I * ratio_K
    I0_a_scaled = ive(0, q1 * a)
    term_sub = (1.0 / I0_a_scaled) * np.exp(-q1 * a) / D_scaled if (I0_a_scaled != 0 and D_scaled != 0) else 0
    return dT_adi_bar_s * (1 - term_sub)

def calc_OLD(tempos, params, T_ini=25, a=0.45):
    res = np.zeros(len(tempos))
    for j, t in enumerate(tempos):
        if t <= 0: res[j] = T_ini; continue
        ln2_t = np.log(2) / t
        soma = sum(V_STEHFEST[i-1] * get_theta_bar_centro_OLD(i * ln2_t, params, a) for i in range(1, 11))
        res[j] = T_ini + soma * ln2_t
    return res

# --- VERSÃO NOVA (COM LEGENDRE 150) ---
def get_theta_bar_centro_NEW(s, params, a):
    dT1, dT2, t1, b1, t2, b2, k_rel, beta_alpha1, beta_alpha2 = params
    alpha1 = beta_alpha1 * 1e-4
    alpha2 = beta_alpha2 * 1e-4
    
    s_arr = np.atleast_1d(s)
    shape = s_arr.shape
    s_f = s_arr.flatten()
    
    # t = -ln(u)/s
    t_n = -np.log(np.maximum(NODES_U.reshape(-1, 1), 1e-15)) / s_f.reshape(1, -1)
    T_v = T_adi_hill_stable(t_n, dT1, dT2, t1, b1, t2, b2)
    dT_adi_bar_f = np.sum(WEIGHTS_U.reshape(-1, 1) * T_v, axis=0) / s_f
    dT_adi_bar = dT_adi_bar_f.reshape(shape)

    q1 = np.sqrt(s / alpha1)
    q2 = np.sqrt(s / alpha2)
    ratio_I = ive(1, q1 * a) / ive(0, q1 * a)
    ratio_K = kve(0, q2 * a) / kve(1, q2 * a)
    flux_ratio = k_rel * np.sqrt(alpha2 / alpha1)
    D_scaled = 1 + flux_ratio * ratio_I * ratio_K
    I0_a_scaled = ive(0, q1 * a)
    
    mask = (I0_a_scaled != 0) & (D_scaled != 0)
    term_sub = np.zeros_like(s)
    term_sub[mask] = (1.0 / I0_a_scaled[mask]) * np.exp(-q1[mask] * a) / D_scaled[mask]
    return dT_adi_bar * (1 - term_sub)

def calc_NEW(tempos, params, T_ini=25, a=0.45):
    tempos = np.atleast_1d(tempos)
    res = np.full(tempos.shape, T_ini, dtype=float)
    mask = tempos > 1e-9
    t_v = tempos[mask]
    if len(t_v) == 0: return res
    ln2_t = np.log(2) / t_v
    S = np.arange(1, 11).reshape(10, 1) * ln2_t.reshape(1, -1)
    theta_bar = get_theta_bar_centro_NEW(S, params, a)
    soma = np.sum(V_STEHFEST.reshape(10, 1) * theta_bar, axis=0)
    res[mask] = T_ini + soma * ln2_t
    return res

if __name__ == "__main__":
    params = [45.0, 40.0, 10.0, 3.0, 25.0, 1.5, 2.9, 40.0, 30.0]
    tempos = np.linspace(0.1, 100, 100)
    
    print("Benchmark: OLD (quad) vs NEW (Legendre 150)")
    s1 = time.time(); r_old = calc_OLD(tempos, params); e_old = time.time() - s1
    s2 = time.time(); r_new = calc_NEW(tempos, params); e_new = time.time() - s2
    
    print(f"Velocidade: OLD={e_old:.3f}s, NEW={e_new:.3f}s (Speedup: {e_old/e_new:.1f}x)")
    print(f"Erro Máximo: {np.max(np.abs(r_old - r_new)):.2e} °C")
    
    plt.plot(tempos, r_old, 'ro', label='OLD (quad)', alpha=0.3)
    plt.plot(tempos, r_new, 'b-', label='NEW (Legendre)')
    plt.legend(); plt.savefig("benchmark_stable.png"); print("Gráfico salvo.")
