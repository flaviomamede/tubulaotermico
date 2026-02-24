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

# Para Gauss-Laguerre (60 pontos)
NODES_GL, WEIGHTS_GL = np.polynomial.laguerre.laggauss(60)

def T_adi_hill(t, dT1, dT2, t1, b1, t2, b2):
    t_safe = np.where(t > 0, t, 1e-9)
    term1 = dT1 * (t_safe**b1) / (t1**b1 + t_safe**b1)
    term2 = dT2 * (t_safe**b2) / (t2**b2 + t_safe**b2)
    return np.where(t > 0, term1 + term2, 0.0)

# --- VERSÃO ANTIGA (BACKEND/MAIN.PY) ---
def get_theta_bar_centro_OLD(s, params, a):
    dT1, dT2, t1, b1, t2, b2, k_rel, beta_alpha1, beta_alpha2 = params
    alpha1 = beta_alpha1 * 1e-4
    alpha2 = beta_alpha2 * 1e-4

    lim_sup = min(2000 / s, 1e5)
    # GARGALO: Integração adaptativa lenta
    dT_adi_bar_s, _ = integrate.quad(
        lambda t: (dT1 * (t**b1) / (t1**b1 + t**b1) + dT2 * (t**b2) / (t2**b2 + t**b2)) * np.exp(-s * t) if t > 0 else 0,
        0, lim_sup, limit=100
    )

    q1 = np.sqrt(s / alpha1)
    q2 = np.sqrt(s / alpha2)
    ratio_I = ive(1, q1 * a) / ive(0, q1 * a)
    ratio_K = kve(0, q2 * a) / kve(1, q2 * a)
    flux_ratio = k_rel * np.sqrt(alpha2 / alpha1)
    D_scaled = 1 + flux_ratio * ratio_I * ratio_K
    I0_a_scaled = ive(0, q1 * a)
    
    if np.isinf(I0_a_scaled) or np.isnan(I0_a_scaled) or D_scaled == 0:
        term_sub = 0
    else:
        term_sub = (1.0 / I0_a_scaled) * np.exp(-q1 * a) / D_scaled
    return dT_adi_bar_s * (1 - term_sub)

def calc_OLD(tempos, params, T_ini=25, a=0.45):
    res = np.zeros(len(tempos))
    # GARGALO: Loop Python sobre os tempos
    for j, t in enumerate(tempos):
        if t <= 0:
            res[j] = T_ini
            continue
        ln2_t = np.log(2) / t
        soma = 0
        for i in range(1, 11):
            soma += V_STEHFEST[i-1] * get_theta_bar_centro_OLD(i * ln2_t, params, a)
        res[j] = T_ini + soma * ln2_t
    return res

# --- VERSÃO NOVA (OTIMIZADA) ---
def calc_NEW(tempos, params, T_ini=25, a=0.45):
    tempos = np.atleast_1d(tempos)
    res = np.full(tempos.shape, T_ini, dtype=float)
    mask = tempos > 0
    t_val = tempos[mask]
    if len(t_val) == 0: return res

    dT1, dT2, t1, b1, t2, b2, k_rel, beta_alpha1, beta_alpha2 = params
    alpha1 = beta_alpha1 * 1e-4
    alpha2 = beta_alpha2 * 1e-4

    # ln(2)/t para Stehfest (Vetorizado sobre o tempo)
    ln2_t = np.log(2) / t_val  # (N_t,)
    
    # s = i * ln(2)/t (Matriz N_stehfest x N_t)
    i_indices = np.arange(1, 11).reshape(10, 1)
    S = i_indices * ln2_t.reshape(1, -1) # (10, N_t)

    # 1. Integração Gauss-Laguerre (Vetorizada sobre s)
    # integral = sum( w_j * T_adi(x_j / s) / s )
    n_nodes = len(NODES_GL)
    nodes_s = NODES_GL.reshape(n_nodes, 1, 1) / S.reshape(1, 10, -1)
    T_nodes = T_adi_hill(nodes_s, dT1, dT2, t1, b1, t2, b2) # (n_nodes, 10, N_t)
    # weight_s: (n_nodes, 1, 1)
    dT_adi_bar = np.sum(WEIGHTS_GL.reshape(n_nodes, 1, 1) * T_nodes / S.reshape(1, 10, -1), axis=0) # (10, N_t)

    # 2. Termo de Difusão (Bessel Vectorized)
    q1 = np.sqrt(S / alpha1)
    q2 = np.sqrt(S / alpha2)
    
    # Ratios de Bessel (ive e kve já são vetorizados pelo numpy)
    ratio_I = ive(1, q1 * a) / ive(0, q1 * a)
    ratio_K = kve(0, q2 * a) / kve(1, q2 * a)
    
    flux_ratio = k_rel * np.sqrt(alpha2 / alpha1)
    D_scaled = 1 + flux_ratio * ratio_I * ratio_K
    I0_a_scaled = ive(0, q1 * a)
    
    # Máscara para evitar divisões por zero ou inf
    mask_valid = (I0_a_scaled != 0) & (D_scaled != 0) & (~np.isinf(I0_a_scaled))
    term_sub = np.zeros_like(S)
    term_sub[mask_valid] = (1.0 / I0_a_scaled[mask_valid]) * np.exp(-q1[mask_valid] * a) / D_scaled[mask_valid]
    
    theta_bar = dT_adi_bar * (1 - term_sub) # (10, N_t)
    
    # 3. Soma de Stehfest Final
    soma = np.sum(V_STEHFEST.reshape(10, 1) * theta_bar, axis=0) # (N_t,)
    res[mask] = T_ini + soma * ln2_t
    return res

# --- BENCHMARK ---
if __name__ == "__main__":
    # Parâmetros de teste (Padrão 9D)
    params = [45.0, 40.0, 10.0, 3.0, 25.0, 1.5, 2.9, 40.0, 30.0]
    # 300 pontos de simulação (mais pesado)
    tempos = np.linspace(0.1, 100, 300)

    print("--- INICIANDO BENCHMARK ---")
    print(f"Testando com {len(tempos)} pontos de tempo...")

    # Versão Antiga
    start = time.time()
    res_old = calc_OLD(tempos, params)
    end_old = time.time() - start
    print(f"Tempo VERSÃO ANTIGA: {end_old:.4f} s")

    # Versão Nova
    start = time.time()
    res_new = calc_NEW(tempos, params)
    end_new = time.time() - start
    print(f"Tempo VERSÃO NOVA:   {end_new:.4f} s")

    # Resultados
    speedup = end_old / end_new
    max_diff = np.max(np.abs(res_old - res_new))
    rmse = np.sqrt(np.mean((res_old - res_new)**2))

    print("\n--- PERFORMANCE ---")
    print(f"Aceleração (Speedup): {speedup:.1f}x")
    
    print("\n--- PRECISÃO ---")
    print(f"Diferença Máxima: {max_diff:.2e} °C")
    print(f"Erro Médio (RMSE): {rmse:.2e} °C")

    if max_diff < 1e-4:
        print("\n[SUCESSO] Os resultados são matematicamente equivalentes!")
    else:
        print("\n[AVISO] Verificou-se uma diferença pequena (> 1e-4), cheque as curvas.")

    # Opcional: Gerar gráfico de comparação
    plt.figure(figsize=(10, 5))
    plt.plot(tempos, res_old, 'o', label='Antiga (quad)', alpha=0.5)
    plt.plot(tempos, res_new, '-', label='Nova (Gauss-Laguerre)', lw=2)
    plt.title(f"Comparação de Modelos (Speedup: {speedup:.1f}x)")
    plt.xlabel("Tempo (h)")
    plt.ylabel("Temperatura (°C)")
    plt.legend()
    plt.grid(True)
    plt.savefig("benchmark_plot.png")
    print("\nGráfico salvo como 'benchmark_plot.png'")
