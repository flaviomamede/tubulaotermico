#!/usr/bin/env python3
"""Diagnóstico completo da regressão - roda direto sem servidor web."""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(__file__)))

import numpy as np
from backend.main import (
    run_otimizacao, calc_temperatura_centro, get_stehfest_V,
    get_theta_bar_centro, _to_beta_scale, ALPHA_SCALE,
    DEFAULT_CHUTE, DEFAULT_BOUNDS_INF, DEFAULT_BOUNDS_SUP
)
from scipy.optimize import least_squares
import scipy.stats as stats

# Usa o mesmo CSV que o usuário — precisa do caminho
# Para teste, usa dados sintéticos de um tubulão 1.4m, T_ini=20
print("="*60)
print("DIAGNÓSTICO DE REGRESSÃO")
print("="*60)

# Parametros "verdadeiros" para gerar dados sintéticos
p_true = [50.0, 45.0, 8.0, 2.5, 40.0, 1.8, 3.5, 35.0, 25.0]
T_ini = 20.0
a = 0.7  # 1.4m diameter / 2

tempos_synth = np.concatenate([
    np.linspace(0.5, 10, 20),   # aquecimento
    np.linspace(12, 40, 15),    # pico
    np.linspace(45, 120, 15),   # resfriamento
])
T_synth = calc_temperatura_centro(tempos_synth, p_true, T_ini=T_ini, a=a)
# Add 0.5°C noise
np.random.seed(42)
T_synth += np.random.normal(0, 0.5, len(T_synth))

print(f"\nDados sintéticos: {len(tempos_synth)} pontos, T range [{T_synth.min():.1f}, {T_synth.max():.1f}]")
print(f"Parâmetros verdadeiros: {p_true}")
print(f"  (alpha1_real = {p_true[7]*ALPHA_SCALE:.4f}, alpha2_real = {p_true[8]*ALPHA_SCALE:.4f})")

# Chute (em escala beta, já convertido)
chute = list(DEFAULT_CHUTE)
bounds_inf = list(DEFAULT_BOUNDS_INF)
bounds_sup = list(DEFAULT_BOUNDS_SUP)

print(f"\nChute: {chute}")
print(f"Bounds inf: {bounds_inf}")
print(f"Bounds sup: {bounds_sup}")

# Rodar otimização direta (sem 2 passos para simplificar o diagnóstico)
def residuals(p):
    return calc_temperatura_centro(tempos_synth, p, T_ini=T_ini, a=a) - T_synth

print("\n--- Rodando least_squares (9 params, todos livres) ---")
res = least_squares(residuals, chute, bounds=(bounds_inf, bounds_sup),
                    method="trf", x_scale="jac", diff_step=1e-6, verbose=2)

print(f"\nStatus: {res.status} — {res.message}")
print(f"nfev: {res.nfev}, njev: {res.njev}")
print(f"Cost (0.5 * sum(r²)): {res.cost:.6f}")
print(f"Optimality: {res.optimality:.2e}")
print(f"\nParâmetros ótimos:")
nomes = ['dT1', 'dT2', 'tau1', 'b1', 'tau2', 'b2', 'k_rel', 'a1(β)', 'a2(β)']
for i, (n, p, pt) in enumerate(zip(nomes, res.x, p_true)):
    at_lb = "⚠LB" if abs(res.x[i] - bounds_inf[i]) < 1e-6 else ""
    at_ub = "⚠UB" if abs(res.x[i] - bounds_sup[i]) < 1e-6 else ""
    print(f"  {n:>8s}: {p:10.4f}  (verdadeiro: {pt:8.4f})  {at_lb}{at_ub}")

# Diagnóstico do Jacobiano
Jac = res.jac
print(f"\nJacobiano shape: {Jac.shape}")
col_norms = np.linalg.norm(Jac, axis=0)
print(f"Normas das colunas: {[f'{v:.4f}' for v in col_norms]}")

# Singular values
U, S, Vt = np.linalg.svd(Jac, full_matrices=False)
print(f"\nValores singulares: {[f'{v:.4e}' for v in S]}")
print(f"Número de condição: {S[0]/S[-1]:.2e}")

# Covariance
n_obs = len(tempos_synth)
p_par = len(res.x)
df = n_obs - p_par
s2 = np.sum(res.fun**2) / df
print(f"\ns² (variância residual): {s2:.6e}")
print(f"RMSE: {np.sqrt(s2):.4f} °C")

FtF = Jac.T @ Jac
print(f"\nDiag(JᵀJ): {[f'{v:.4f}' for v in np.diag(FtF)]}")
print(f"cond(JᵀJ): {np.linalg.cond(FtF):.2e}")

try:
    FtF_inv = np.linalg.inv(FtF)
    Sigma = s2 * FtF_inv
    SE = np.sqrt(np.clip(np.diag(Sigma), 0, None))
    print(f"\nSE (inv direto): {[f'{v:.4f}' for v in SE]}")
except np.linalg.LinAlgError:
    print("\n⚠ JᵀJ é singular, usando pinv")

FtF_inv_p = np.linalg.pinv(FtF, rcond=1e-10)
Sigma_p = s2 * FtF_inv_p
SE_p = np.sqrt(np.clip(np.diag(Sigma_p), 0, None))
print(f"SE (pinv): {[f'{v:.4f}' for v in SE_p]}")

# Agora testar com o preconditioning do código original
scale_factors = col_norms.copy()
scale_factors[scale_factors < 1e-12] = 1.0
F_scaled = Jac / scale_factors
FtF_s = F_scaled.T @ F_scaled
print(f"\ncond(F_scaled.T @ F_scaled): {np.linalg.cond(FtF_s):.2e}")
FtF_s_inv = np.linalg.pinv(FtF_s, rcond=1e-5)
D_inv = np.diag(1.0 / scale_factors)
FtF_inv_pc = D_inv @ FtF_s_inv @ D_inv
Sigma_pc = s2 * FtF_inv_pc
SE_pc = np.sqrt(np.clip(np.diag(Sigma_pc), 0, None))
print(f"SE (preconditioned pinv): {[f'{v:.4f}' for v in SE_pc]}")

print(f"\n{'='*60}")
print("Resumo: se os SEs acima são razoáveis (< 10), a estatística funciona.")
print("Se explodem (> 1e6), há multicolinearidade entre os parâmetros.")
print(f"{'='*60}")
