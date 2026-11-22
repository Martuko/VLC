# -*- coding: utf-8 -*-
"""
run_all.py
==========
Ejecuta todos los experimentos del proyecto "Simulación de un sistema VLC en un
entorno industrial hostil" y guarda figuras (PNG) y resultados (CSV) en ./resultados

Experimentos
------------
E1  Plano de la nave e iluminancia (verifica que la instalación sirve para iluminar)
E2  Caso de referencia: potencia recibida, SNR y dispersión de retardo RMS
E3  Mapas de SNR para 6 escenarios (ideal, polvo, temperatura, luz ambiente,
    obstáculos y combinado) + cobertura confiable (BER <= 1e-6)
E4  Polvo en suspensión: BER vs visibilidad para varias tasas de bits
E5  Suciedad en las ópticas: cobertura vs días sin limpieza
E6  Temperatura: SNR vs temperatura ambiente (sensibilidad del LED y del receptor)
E7  Tasa de bits e ISI: cobertura vs Rb
E8  Obstáculos móviles (Monte Carlo): disponibilidad vs número de operarios,
    para distintas configuraciones (FOV, número de luminarias) y Wi-Fi / híbrido
E9  Comparación VLC vs Wi-Fi 2,4 GHz: cobertura vs EMI y vs polvo
E10 Diseño mejorado: escenario combinado antes y después de las recomendaciones

Uso:  python run_all.py            (≈ 2-4 min en un PC normal)
"""
import os, time, csv
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.colors import ListedColormap, BoundaryNorm
from vlc_core import *

OUT = "resultados"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({"font.size": 10, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 110, "savefig.dpi": 160, "savefig.bbox": "tight"})
RNG = np.random.default_rng(2023)
BER_TH = 1e-6
t_start = time.time()

def log(msg):
    print(f"[{time.time() - t_start:6.1f} s] {msg}", flush=True)

def beta_of(V):
    return 0.0 if np.isinf(V) else 3.912 / V

def coverage(ber):
    return 100.0 * np.mean(ber <= BER_TH)

def draw_layout(ax, p, boxes_dyn=None, show_leds=True, show_ap=False, show_static=True):
    for name, b in (STATIC_OBSTACLES.items() if show_static else []):
        ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2], fc="0.55", ec="k", lw=0.8, alpha=0.85))
    if boxes_dyn is not None:
        for b in boxes_dyn:
            c = "tab:orange" if (b[1] - b[0]) < 0.6 else "tab:red"
            ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2], fc=c, ec="k", lw=0.5))
    if show_leds:
        L = p.leds
        ax.plot(L[:, 0], L[:, 1], "*", ms=11, mfc="gold", mec="k", label="Luminaria VLC")
    if show_ap:
        w = WifiParams()
        ax.plot(w.ap_pos[0], w.ap_pos[1], "^", ms=9, mfc="tab:cyan", mec="k", label="AP Wi-Fi")
    ax.set_xlim(0, p.L); ax.set_ylim(0, p.W); ax.set_aspect("equal")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.grid(False)

csv_rows = []   # resumen general

# ===========================================================================
# Geometría base
# ===========================================================================
P0 = Params()
SB = static_boxes()
rx_map, xs, ys = receiver_grid(P0, 0.25)
free_map = ~inside_any(rx_map[:, :2], SB)          # puntos fuera de las máquinas
nx, ny = len(xs), len(ys)
log(f"Grilla de mapas: {len(rx_map)} puntos ({free_map.sum()} en área libre)")

# ===========================================================================
# E1  Plano + iluminancia
# ===========================================================================
E = illuminance(P0, rx_map)
fig, axs = plt.subplots(1, 2, figsize=(11, 4.6))
draw_layout(axs[0], P0, show_ap=True)
for name, b in STATIC_OBSTACLES.items():
    axs[0].text((b[0] + b[1]) / 2, b[3] - 0.15, f"{name} ({b[4]:.1f} m)", ha="center", va="top", fontsize=6.5)
axs[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.13), ncol=2, fontsize=8)
axs[0].set_title("Planta de la nave (15 × 15 × 6 m)")
im = axs[1].imshow(E.reshape(ny, nx), origin="lower", extent=[0, P0.L, 0, P0.W], cmap="inferno")
cs = axs[1].contour(xs, ys, E.reshape(ny, nx), levels=[200, 300, 400], colors="w", linewidths=0.8)
axs[1].clabel(cs, fmt="%d lx", fontsize=7)
plt.colorbar(im, ax=axs[1], label="Iluminancia [lx]")
axs[1].set_title(f"Iluminancia directa en el plano de 1 m\nmedia {E.mean():.0f} lx · U0 = {E.min()/E.mean():.2f}")
axs[1].set_xlabel("x [m]"); axs[1].set_ylabel("y [m]"); axs[1].grid(False)
fig.savefig(f"{OUT}/E1_planta_iluminancia.png"); plt.close(fig)
csv_rows.append(["E1", "Iluminancia media [lx]", f"{E.mean():.0f}"])
csv_rows.append(["E1", "Uniformidad U0 = Emin/Emed", f"{E.min()/E.mean():.2f}"])
log("E1 listo")

# ===========================================================================
# Canales base (mapas): sin obstáculos y con obstáculos fijos, aire limpio y polvo V = 20 m
# ===========================================================================
V_DUST = 20.0
betas_map = (0.0, beta_of(V_DUST))
h0, d0 = los_channel(P0, rx_map)
hist0, tb = reflection_histogram(P0, rx_map, None, betas_map)
log("Canal sin obstáculos listo")
hS, dS = los_channel(P0, rx_map, SB)
histS, _ = reflection_histogram(P0, rx_map, SB, betas_map)
log("Canal con obstáculos fijos listo")

# ===========================================================================
# E2  Caso de referencia
# ===========================================================================
r0 = performance(P0, h0, d0, hist0[0.0], tb)
tau = rms_delay_spread(h0, d0, hist0[0.0], tb)
fig, axs = plt.subplots(1, 3, figsize=(15, 4.3))
for ax, data, lab, cm in ((axs[0], 10 * np.log10(r0["P_dc"] * 1e3), "Potencia óptica recibida [dBm]", "viridis"),
                          (axs[1], r0["snr_db"], "SNR [dB]", "viridis"),
                          (axs[2], tau * 1e9, "Dispersión de retardo RMS [ns]", "magma")):
    im = ax.imshow(data.reshape(ny, nx), origin="lower", extent=[0, P0.L, 0, P0.W], cmap=cm)
    plt.colorbar(im, ax=ax, label=lab)
    ax.plot(P0.leds[:, 0], P0.leds[:, 1], "*", ms=9, mfc="gold", mec="k")
    ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]"); ax.grid(False)
axs[0].set_title(f"Potencia recibida ({10*np.log10(r0['P_dc'].min()*1e3):.1f} a {10*np.log10(r0['P_dc'].max()*1e3):.1f} dBm)")
axs[1].set_title(f"SNR a 50 Mb/s ({r0['snr_db'].min():.1f} a {r0['snr_db'].max():.1f} dB)")
axs[2].set_title(f"τ_RMS (máx. {tau.max()*1e9:.1f} ns)")
fig.suptitle("E2 · Caso de referencia: nave sin obstáculos, aire limpio, 25 °C", y=1.02)
fig.savefig(f"{OUT}/E2_referencia_potencia_snr_retardo.png"); plt.close(fig)
csv_rows += [["E2", "Potencia recibida min/max [dBm]", f"{10*np.log10(r0['P_dc'].min()*1e3):.1f} / {10*np.log10(r0['P_dc'].max()*1e3):.1f}"],
             ["E2", "SNR min/mediana/max [dB]", f"{r0['snr_db'].min():.1f} / {np.median(r0['snr_db']):.1f} / {r0['snr_db'].max():.1f}"],
             ["E2", "tau_RMS max [ns]", f"{tau.max()*1e9:.1f}"]]

# Respuesta impulsiva en dos puntos (centro y esquina)
fig, ax = plt.subplots(figsize=(7, 3.6))
for (x, y, c) in ((7.5, 7.5, "tab:blue"), (0.625, 0.625, "tab:red")):
    i = np.argmin((rx_map[:, 0] - x) ** 2 + (rx_map[:, 1] - y) ** 2)
    hr = hist0[0.0][i] / (tb[1] - tb[0])
    ax.plot(tb * 1e9, hr, c, lw=1, label=f"Reflexiones en ({rx_map[i,0]:.1f}, {rx_map[i,1]:.1f}) m")
    for j in range(len(P0.leds)):
        if h0[i, j] > 0:
            ax.vlines(d0[i, j] / C_LUZ * 1e9, 0, hr.max() * 1.15, colors=c, lw=1.5, alpha=0.6)
ax.set_xlim(10, 90); ax.set_xlabel("Tiempo [ns]"); ax.set_ylabel("h(t) [1/s]")
ax.set_title("Respuesta impulsiva (líneas verticales = llegadas LOS de cada luminaria)")
ax.legend(fontsize=8)
fig.savefig(f"{OUT}/E2_respuesta_impulsiva.png"); plt.close(fig)
log("E2 listo")

# ===========================================================================
# E3  Seis escenarios
# ===========================================================================
dyn_map = random_dynamic_boxes(RNG, 12, 3, P0)
hD, dD = los_channel(P0, rx_map, np.vstack([SB, dyn_map]))   # LOS con obstáculos fijos+móviles
scen = {
    "S0 · Ideal":                (P0, h0, d0, hist0[0.0]),
    "S1 · Polvo (V=20 m, 30 días)": (replace(P0, visibility=V_DUST, days_dirty=30), h0, d0, hist0[beta_of(V_DUST)]),
    "S2 · Temperatura 70 °C":    (replace(P0, T_amb=70), h0, d0, hist0[0.0]),
    "S3 · Luz ambiente intensa": (replace(P0, I_bg=5100e-6), h0, d0, hist0[0.0]),
    "S4 · Obstáculos":           (P0, hD, dD, histS[0.0]),
    "S5 · Combinado":            (replace(P0, visibility=V_DUST, days_dirty=30, T_amb=70, I_bg=5100e-6),
                                  hD, dD, histS[beta_of(V_DUST)]),
}
res_scen = {}
fig, axs = plt.subplots(2, 3, figsize=(15.5, 9.6))
cmap = plt.get_cmap("RdYlGn")
for ax, (name, (p, h, d, hi)) in zip(axs.ravel(), scen.items()):
    r = performance(p, h, d, hi, tb)
    with_obs = ("Obst" in name) or ("Comb" in name)
    area = free_map if with_obs else np.ones(len(free_map), bool)
    snr = np.where(area, r["snr_db"], np.nan)
    cov = coverage(r["ber"][area])
    res_scen[name] = (r, cov, area)
    im = ax.imshow(np.clip(snr, -10, 50).reshape(ny, nx), origin="lower", extent=[0, P0.L, 0, P0.W],
                   cmap=cmap, vmin=-10, vmax=50)
    ax.contour(xs, ys, np.nan_to_num(snr, nan=99).reshape(ny, nx), levels=[SNR_TH_OOK_DB], colors="k", linewidths=1.2)
    draw_layout(ax, P0, dyn_map if with_obs else None, show_static=with_obs)
    ax.set_title(f"{name}\ncobertura {cov:.1f} % · SNR med. {np.nanmedian(snr):.1f} dB", fontsize=9.5)
    csv_rows.append(["E3", f"{name} cobertura [%]", f"{cov:.1f}"])
    csv_rows.append(["E3", f"{name} SNR mediana [dB]", f"{np.nanmedian(snr):.1f}"])
fig.colorbar(im, ax=axs, shrink=0.8, label="SNR [dB]  (línea negra = 13,5 dB → BER = 10⁻⁶)")
fig.suptitle("E3 · SNR a 50 Mb/s en seis escenarios · cobertura = % del área con BER ≤ 10⁻⁶\n(S4–S5: gris = máquinas, naranjo = operarios, rojo = grúas horquilla)")
fig.savefig(f"{OUT}/E3_mapas_snr_escenarios.png"); plt.close(fig)

# CDF de BER por escenario
fig, ax = plt.subplots(figsize=(7.5, 4.2))
for name, (r, cov, area) in res_scen.items():
    b = np.sort(np.log10(np.maximum(r["ber"][area], 1e-30)))
    ax.plot(b, np.arange(1, len(b) + 1) / len(b) * 100, lw=1.6, label=name)
ax.axvline(-6, color="k", ls="--", lw=1); ax.text(-5.9, 5, "BER = 10⁻⁶", fontsize=8)
ax.set_xlim(-30, 0); ax.set_xlabel("log10(BER)"); ax.set_ylabel("% del área con BER ≤ x")
ax.set_title("E3 · Distribución acumulada de la BER en el área de trabajo"); ax.legend(fontsize=8)
fig.savefig(f"{OUT}/E3_cdf_ber.png"); plt.close(fig)
log("E3 listo")

# ===========================================================================
# Grilla gruesa para barridos (0,5 m) con obstáculos fijos
# ===========================================================================
rx_s, _, _ = receiver_grid(P0, 0.5)
rx_s = rx_s[~inside_any(rx_s[:, :2], SB)]
V_list = [np.inf, 200, 100, 50, 30, 20, 15, 10, 7, 5]
betas_s = tuple(sorted(set(beta_of(v) for v in V_list)))
h_s, d_s = los_channel(P0, rx_s, SB)
hist_s, tb = reflection_histogram(P0, rx_s, SB, betas_s)
h_s0, d_s0 = los_channel(P0, rx_s)
hist_s0, _ = reflection_histogram(P0, rx_s, None, betas_s)
log(f"Grilla de barridos lista ({len(rx_s)} puntos)")

# ===========================================================================
# E4  Polvo: BER vs visibilidad para varias tasas
# ===========================================================================
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
Vplot = [v for v in V_list if not np.isinf(v)]
rows_e4 = []
for Rb, c in ((10e6, "tab:green"), (50e6, "tab:blue"), (100e6, "tab:red")):
    med, p95, cov = [], [], []
    for v in Vplot:
        p = replace(P0, visibility=v, Rb=Rb)
        r = performance(p, h_s0, d_s0, hist_s0[beta_of(v)], tb)
        lb = np.log10(np.maximum(r["ber"], 1e-40))
        med.append(np.median(lb)); p95.append(np.percentile(lb, 95)); cov.append(coverage(r["ber"]))
        rows_e4.append([Rb / 1e6, v, np.median(r["snr_db"]), np.percentile(r["snr_db"], 5), coverage(r["ber"])])
    axs[0].semilogx(Vplot, med, "o-", c=c, label=f"{Rb/1e6:.0f} Mb/s · mediana")
    axs[0].semilogx(Vplot, p95, "s--", c=c, alpha=0.7, label=f"{Rb/1e6:.0f} Mb/s · peor 5 %")
    axs[1].semilogx(Vplot, cov, "o-", c=c, label=f"{Rb/1e6:.0f} Mb/s")
axs[0].axhline(-6, color="k", ls=":", lw=1)
axs[0].set_ylim(-30, 0); axs[0].set_xlabel("Visibilidad V [m] (menor = más polvo)"); axs[0].set_ylabel("log10(BER)")
axs[0].set_title("E4 · BER vs concentración de polvo (sin obstáculos)"); axs[0].legend(fontsize=7, ncol=2)
axs[0].invert_xaxis(); axs[1].invert_xaxis()
axs[1].set_xlabel("Visibilidad V [m]"); axs[1].set_ylabel("Cobertura con BER ≤ 10⁻⁶ [%]")
axs[1].set_title("E4 · Cobertura vs polvo"); axs[1].legend(fontsize=8)
for ax in axs:
    ax.set_xticks([200, 100, 50, 20, 10, 5]); ax.set_xticklabels(["200", "100", "50", "20", "10", "5"])
fig.savefig(f"{OUT}/E4_polvo_ber.png"); plt.close(fig)
with open(f"{OUT}/E4_polvo.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["Rb_Mbps", "visibilidad_m", "SNR_mediana_dB", "SNR_p5_dB", "cobertura_%"]); wr.writerows(rows_e4)
log("E4 listo")

# ===========================================================================
# E5  Suciedad en ópticas: cobertura vs días sin limpieza
# ===========================================================================
days = np.arange(0, 121, 5)
fig, ax = plt.subplots(figsize=(7.5, 4.2))
rows_e5 = []
for v, c in ((np.inf, "tab:green"), (50, "tab:blue"), (20, "tab:orange"), (10, "tab:red")):
    cov = []
    for dd in days:
        p = replace(P0, visibility=v, days_dirty=dd)
        r = performance(p, h_s, d_s, hist_s[beta_of(v)], tb)
        cov.append(coverage(r["ber"])); rows_e5.append([v, dd, coverage(r["ber"]), np.median(r["snr_db"])])
    lab = "aire limpio" if np.isinf(v) else f"polvo V = {v:.0f} m"
    ax.plot(days, cov, "o-", ms=3, c=c, label=lab)
ax.set_xlabel("Días sin limpiar luminarias y receptores"); ax.set_ylabel("Cobertura con BER ≤ 10⁻⁶ [%]")
ax.set_title("E5 · Efecto de la suciedad acumulada (con máquinas, 50 Mb/s)"); ax.legend(fontsize=8)
fig.savefig(f"{OUT}/E5_suciedad_limpieza.png"); plt.close(fig)
with open(f"{OUT}/E5_suciedad.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["visibilidad_m", "dias", "cobertura_%", "SNR_mediana_dB"]); wr.writerows(rows_e5)
log("E5 listo")

# ===========================================================================
# E6  Temperatura
# ===========================================================================
Ts = np.arange(0, 91, 5)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
ref = performance(P0, h_s0, d_s0, hist_s0[0.0], tb)
for kT, c, ls in ((0.002, "tab:green", "-"), (0.004, "tab:blue", "-"), (0.010, "tab:red", "-")):
    dsnr = [np.median(performance(replace(P0, T_amb=T, k_T_led=kT), h_s0, d_s0, hist_s0[0.0], tb)["snr_db"])
            - np.median(ref["snr_db"]) for T in Ts]
    axs[0].plot(Ts, dsnr, ls, c=c, label=f"LED {kT*100:.1f} %/°C")
dsnr = [np.median(performance(replace(P0, T_amb=T, T_derate_C=60), h_s0, d_s0, hist_s0[0.0], tb)["snr_db"])
        - np.median(ref["snr_db"]) for T in Ts]
axs[0].plot(Ts, dsnr, "--", c="k", label="0,4 %/°C + protección térmica >60 °C")
axs[0].set_xlabel("Temperatura ambiente [°C]"); axs[0].set_ylabel("Variación de la SNR mediana [dB]")
axs[0].set_title("E6 · Efecto de la temperatura sobre la SNR"); axs[0].legend(fontsize=8)
# descomposición del ruido en un punto típico (centro de la nave) a 25 y 85 °C
i = np.argmin((rx_s[:, 0] - 7.25) ** 2 + (rx_s[:, 1] - 6.25) ** 2)
labels, sh, th, dk = [], [], [], []
for T in (25, 50, 70, 85):
    p = replace(P0, T_amb=T)
    r = performance(p, h_s0[i:i+1], d_s0[i:i+1], hist_s0[0.0][i:i+1], tb)
    I_dark = p.I_dark0 * 2 ** ((T - 25) / 10)
    s_dark = 2 * Q_E * I_dark * p.Rb
    labels.append(f"{T} °C"); sh.append(float(r["s_shot"][0]) - s_dark); th.append(float(r["s_th"])); dk.append(s_dark)
x = np.arange(len(labels))
axs[1].bar(x - 0.27, sh, 0.27, label="Shot (señal + luz ambiente)")
axs[1].bar(x, th, 0.27, label="Térmico")
axs[1].bar(x + 0.27, dk, 0.27, label="Corriente oscura")
axs[1].set_yscale("log"); axs[1].set_xticks(x); axs[1].set_xticklabels(labels)
axs[1].set_ylabel("Varianza de ruido [A²]"); axs[1].set_title("E6 · Componentes de ruido vs temperatura (centro, 50 Mb/s)")
axs[1].legend(fontsize=8, loc="upper center", bbox_to_anchor=(0.5, -0.1), ncol=3)
fig.savefig(f"{OUT}/E6_temperatura.png"); plt.close(fig)
csv_rows.append(["E6", "ΔSNR mediana a 70 °C, LED 0,4 %/°C [dB]",
                 f"{np.median(performance(replace(P0, T_amb=70), h_s0, d_s0, hist_s0[0.0], tb)['snr_db']) - np.median(ref['snr_db']):.2f}"])
csv_rows.append(["E6", "ΔSNR mediana a 85 °C con protección térmica [dB]", f"{dsnr[-2]:.2f}"])
log("E6 listo")

# ===========================================================================
# E7  Tasa de bits e ISI
# ===========================================================================
Rbs = np.array([1, 2, 5, 10, 20, 50, 75, 100, 150, 200]) * 1e6
fig, ax = plt.subplots(figsize=(7.5, 4.2))
rows_e7 = []
cases = (("Ideal", P0, h_s0, d_s0, hist_s0[0.0], "tab:green"),
         ("Polvo V=20 m + 30 días", replace(P0, visibility=20, days_dirty=30), h_s0, d_s0, hist_s0[beta_of(20)], "tab:orange"),
         ("Combinado (fijos, 70 °C, sol)", replace(P0, visibility=20, days_dirty=30, T_amb=70, I_bg=5100e-6),
          h_s, d_s, hist_s[beta_of(20)], "tab:red"))
for name, p, h, d, hi, c in cases:
    cov = []
    for Rb in Rbs:
        r = performance(replace(p, Rb=Rb), h, d, hi, tb)
        cov.append(coverage(r["ber"])); rows_e7.append([name, Rb / 1e6, coverage(r["ber"]), np.median(r["snr_db"])])
    ax.semilogx(Rbs / 1e6, cov, "o-", c=c, label=name)
ax.set_xlabel("Tasa de bits OOK [Mb/s]"); ax.set_ylabel("Cobertura con BER ≤ 10⁻⁶ [%]")
ax.set_title("E7 · Cobertura vs tasa de bits (ruido ∝ B y B³, ISI por reflexiones)"); ax.legend(fontsize=8)
fig.savefig(f"{OUT}/E7_tasa_bits.png"); plt.close(fig)
with open(f"{OUT}/E7_tasa_bits.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["escenario", "Rb_Mbps", "cobertura_%", "SNR_mediana_dB"]); wr.writerows(rows_e7)
log("E7 listo")

# ===========================================================================
# E8  Obstáculos móviles (Monte Carlo): VLC (varias configuraciones), Wi-Fi e híbrido
# ===========================================================================
N_MC = 300
rx_mc = rx_s[RNG.choice(len(rx_s), 250, replace=False)]
W0 = WifiParams(emi_dB=20.0)                  # EMI moderada en la nave
ap = np.array(W0.ap_pos)
d_ap = np.linalg.norm(rx_mc - ap, axis=1)
SF = RNG.standard_normal(len(rx_mc))          # sombreamiento fijo por ubicación (normalizado)
N_dbm = wifi_noise_dBm(W0)

def vlc_config(name, **kw):
    p = replace(P0, **kw)
    hist, _ = reflection_histogram(p, rx_mc, SB, (0.0,))
    return name, p, hist[0.0]

configs = [vlc_config("VLC · 9 luminarias, FOV 60°"),
           vlc_config("VLC · 9 luminarias, FOV 80°", fov_deg=80.0),
           vlc_config("VLC · 9 luminarias, FOV 40°", fov_deg=40.0),
           vlc_config("VLC · 4 luminarias (misma iluminancia), FOV 60°",
                      led_xy=((3.75, 3.75), (3.75, 11.25), (11.25, 3.75), (11.25, 11.25)), P_opt=45.0 * 9 / 4)]
n_workers_list = [0, 5, 10, 20, 30, 40]
avail = {c[0]: [] for c in configs}
avail["Wi-Fi 2,4 GHz (EMI +20 dB)"] = []
avail["Híbrido VLC (9 lum., FOV 60°) + Wi-Fi"] = []
for nw in n_workers_list:
    ok = {k: np.zeros(len(rx_mc)) for k in avail}
    for it in range(N_MC):
        dyn = random_dynamic_boxes(RNG, nw, nw // 10, P0)
        boxes = np.vstack([SB, dyn])
        vlc_ok_main = None
        for name, p, hi in configs:
            h, d = los_channel(p, rx_mc, boxes)
            r = performance(p, h, d, hi, tb)
            good = r["ber"] <= BER_TH
            ok[name] += good
            if vlc_ok_main is None:
                vlc_ok_main = good
        los = ~segment_blocked(rx_mc, np.repeat(ap[None, :], len(rx_mc), 0), boxes)
        pl = wifi_pathloss(W0, d_ap, los)
        sf = np.where(los, W0.sigma_los, W0.sigma_nlos) * SF
        fad = wifi_fading_power(RNG, los, 1)[:, 0]
        snr_w = W0.Pt_dBm + 2 * W0.G_ant_dBi - pl - sf + 10 * np.log10(fad) - N_dbm
        wgood = snr_w >= SNR_TH_BPSK_DB
        ok["Wi-Fi 2,4 GHz (EMI +20 dB)"] += wgood
        ok["Híbrido VLC (9 lum., FOV 60°) + Wi-Fi"] += (wgood | vlc_ok_main)
    for k in avail:
        # % de ubicaciones con disponibilidad >= 99 %
        avail[k].append(100 * np.mean(ok[k] / N_MC >= 0.99))
    log(f"E8: {nw} operarios listo")

fig, ax = plt.subplots(figsize=(8.5, 4.6))
styles = ["o-", "s-", "^-", "v-", "D--", "P-"]
colors = ["tab:blue", "tab:green", "tab:purple", "tab:gray", "tab:cyan", "k"]
for (k, v), st, c in zip(avail.items(), styles, colors):
    ax.plot(n_workers_list, v, st, c=c, label=k)
ax.set_xlabel("Operarios en la nave (+1 grúa horquilla por cada 10)")
ax.set_ylabel("% de ubicaciones con disponibilidad ≥ 99 %")
ax.set_title("E8 · Bloqueo por obstáculos móviles (Monte Carlo, 300 realizaciones)")
ax.legend(fontsize=7.5, loc="lower left"); ax.set_ylim(0, 102)
fig.savefig(f"{OUT}/E8_obstaculos_montecarlo.png"); plt.close(fig)
with open(f"{OUT}/E8_obstaculos.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["configuracion"] + [f"{n}_operarios" for n in n_workers_list])
    for k, v in avail.items():
        wr.writerow([k] + [f"{x:.1f}" for x in v])
log("E8 listo")

# ===========================================================================
# E9  Comparación VLC vs Wi-Fi: EMI y polvo
# ===========================================================================
rx_c = rx_s
d_apc = np.linalg.norm(rx_c - ap, axis=1)
los_c = ~segment_blocked(rx_c, np.repeat(ap[None, :], len(rx_c), 0), SB)
SFc = RNG.standard_normal(len(rx_c))
FAD = wifi_fading_power(RNG, los_c, 2000)

def wifi_cov(emi):
    w = replace(W0, emi_dB=emi)
    pl = wifi_pathloss(w, d_apc, los_c)
    sf = np.where(los_c, w.sigma_los, w.sigma_nlos) * SFc
    snr = w.Pt_dBm + 2 * w.G_ant_dBi - pl - sf
    snr_inst = snr[:, None] + 10 * np.log10(FAD) - wifi_noise_dBm(w)
    av = np.mean(snr_inst >= SNR_TH_BPSK_DB, axis=1)
    ber_mean = np.mean(qfunc(np.sqrt(2 * 10 ** (snr_inst / 10))), axis=1)
    return 100 * np.mean(av >= 0.99), np.median(snr - wifi_noise_dBm(w)), np.median(ber_mean)

emis = [0, 10, 20, 30, 40, 50]
wifi_e = [wifi_cov(e) for e in emis]
vlc_ref = performance(P0, h_s, d_s, hist_s[0.0], tb)
vlc_cov_emi = [coverage(vlc_ref["ber"])] * len(emis)
fig, axs = plt.subplots(1, 2, figsize=(12, 4.3))
axs[0].plot(emis, [w[0] for w in wifi_e], "D-", c="tab:cyan", label="Wi-Fi 2,4 GHz (InF-SH, BPSK 20 MHz)")
axs[0].plot(emis, vlc_cov_emi, "o-", c="tab:blue", label="VLC 50 Mb/s (inmune a EMI)")
axs[0].set_xlabel("Aumento del piso de ruido RF por EMI industrial [dB]")
axs[0].set_ylabel("Cobertura confiable [%]"); axs[0].set_ylim(0, 102)
axs[0].set_title("E9a · Interferencia electromagnética (motores, VFD, soldadura)"); axs[0].legend(fontsize=8)
vlc_cov_v, wifi_cov_v = [], []
for v in Vplot:
    r = performance(replace(P0, visibility=v), h_s, d_s, hist_s[beta_of(v)], tb)
    vlc_cov_v.append(coverage(r["ber"]))
wifi_cov_v = [wifi_cov(20)[0]] * len(Vplot)
axs[1].semilogx(Vplot, vlc_cov_v, "o-", c="tab:blue", label="VLC 50 Mb/s")
axs[1].semilogx(Vplot, wifi_cov_v, "D-", c="tab:cyan", label="Wi-Fi (EMI +20 dB; el polvo no la afecta)")
axs[1].invert_xaxis(); axs[1].set_xticks([200, 100, 50, 20, 10, 5]); axs[1].set_xticklabels(["200", "100", "50", "20", "10", "5"])
axs[1].set_xlabel("Visibilidad V [m]"); axs[1].set_ylabel("Cobertura confiable [%]"); axs[1].set_ylim(0, 102)
axs[1].set_title("E9b · Polvo en suspensión"); axs[1].legend(fontsize=8)
fig.savefig(f"{OUT}/E9_comparacion_wifi.png"); plt.close(fig)
with open(f"{OUT}/E9_wifi.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["EMI_dB", "wifi_cobertura_%", "wifi_SNR_mediana_dB", "wifi_BER_media_mediana", "vlc_cobertura_%"])
    for e, w in zip(emis, wifi_e):
        wr.writerow([e, f"{w[0]:.1f}", f"{w[1]:.1f}", f"{w[2]:.2e}", f"{vlc_cov_emi[0]:.1f}"])
csv_rows.append(["E9", "VLC SNR mediana (con máquinas) [dB]", f"{np.median(vlc_ref['snr_db']):.1f}"])
csv_rows.append(["E9", "VLC BER mediana (con máquinas)", f"{np.median(vlc_ref['ber']):.1e}"])
log("E9 listo")

# ===========================================================================
# E10 Diseño mejorado vs base en el escenario combinado (Monte Carlo)
# ===========================================================================
hard = dict(visibility=20.0, T_amb=70.0)
base = replace(P0, days_dirty=30, I_bg=5100e-6, **hard)
improved = replace(P0, days_dirty=7, I_bg=5100e-6 * 0.3, Ts_filter=0.85, Rb=10e6, **hard)
designs = {"Diseño base": base, "Diseño mejorado": improved}
hist_d = {k: reflection_histogram(p, rx_mc, SB, (beta_of(20.0),))[0][beta_of(20.0)] for k, p in designs.items()}
okD = {k: np.zeros(len(rx_mc)) for k in list(designs) + ["Mejorado + respaldo Wi-Fi"]}
snr_store = {k: [] for k in designs}
for it in range(N_MC):
    dyn = random_dynamic_boxes(RNG, 12, 3, P0)
    boxes = np.vstack([SB, dyn])
    for k, p in designs.items():
        h, d = los_channel(p, rx_mc, boxes)
        r = performance(p, h, d, hist_d[k], tb)
        okD[k] += r["ber"] <= BER_TH
        snr_store[k].append(r["snr_db"])
        if k == "Diseño mejorado":
            vgood = r["ber"] <= BER_TH
    los = ~segment_blocked(rx_mc, np.repeat(ap[None, :], len(rx_mc), 0), boxes)
    pl = wifi_pathloss(W0, d_ap, los)
    sf = np.where(los, W0.sigma_los, W0.sigma_nlos) * SF
    fad = wifi_fading_power(RNG, los, 1)[:, 0]
    snr_w = W0.Pt_dBm + 2 * W0.G_ant_dBi - pl - sf + 10 * np.log10(fad) - N_dbm
    okD["Mejorado + respaldo Wi-Fi"] += vgood | (snr_w >= SNR_TH_BPSK_DB)
e10 = {k: 100 * np.mean(v / N_MC >= 0.99) for k, v in okD.items()}
fig, axs = plt.subplots(1, 2, figsize=(12, 4.2))
for k, c in (("Diseño base", "tab:red"), ("Diseño mejorado", "tab:green")):
    s = np.sort(np.concatenate(snr_store[k]))
    axs[0].plot(s, np.arange(1, len(s) + 1) / len(s) * 100, c=c, lw=1.8, label=k)
axs[0].axvline(SNR_TH_OOK_DB, color="k", ls="--", lw=1); axs[0].text(SNR_TH_OOK_DB + 0.5, 50, "13,5 dB\n(BER 10⁻⁶)", fontsize=8)
axs[0].set_xlim(-20, 50); axs[0].set_xlabel("SNR [dB]"); axs[0].set_ylabel("CDF [%]")
axs[0].set_title("E10 · SNR en el escenario combinado"); axs[0].legend(fontsize=8)
ks = list(e10); vals = [e10[k] for k in ks]
bars = axs[1].bar(range(len(ks)), vals, color=["tab:red", "tab:green", "tab:olive"])
for b_, v in zip(bars, vals):
    axs[1].text(b_.get_x() + b_.get_width() / 2, v + 1, f"{v:.0f} %", ha="center")
axs[1].set_xticks(range(len(ks))); axs[1].set_xticklabels(ks, fontsize=8.5)
axs[1].set_ylabel("% de ubicaciones con disponibilidad ≥ 99 %"); axs[1].set_ylim(0, 110)
axs[1].set_title("E10 · Disponibilidad: base vs recomendaciones")
fig.savefig(f"{OUT}/E10_diseno_mejorado.png"); plt.close(fig)
for k, v in e10.items():
    csv_rows.append(["E10", f"{k} · ubicaciones con disponibilidad ≥99 % [%]", f"{v:.1f}"])
log("E10 listo")

with open(f"{OUT}/resumen_resultados.csv", "w", newline="") as f:
    wr = csv.writer(f); wr.writerow(["experimento", "indicador", "valor"]); wr.writerows(csv_rows)
for r_ in csv_rows:
    print(" | ".join(r_))
print("\nE8 disponibilidad:"); [print(k, [f"{x:.0f}" for x in v]) for k, v in avail.items()]
print("\nE9 Wi-Fi (EMI, cobertura, SNR med, BER med):"); [print(e, w) for e, w in zip(emis, wifi_e)]
print("E9 VLC cobertura vs V:", list(zip(Vplot, [f"{c:.1f}" for c in vlc_cov_v])))
log("Terminado")
