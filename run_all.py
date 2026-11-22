# -*- coding: utf-8 -*-
"""
run_all.py
==========
Ejecuta los 10 experimentos (E1–E10) del simulador VLC industrial mejorado en
una nave metalmecánica pesada con alta EMI (horno de inducción, prensa, CNCs
con VFD, puesto de soldadura, polvo en suspensión, suciedad óptica y tráfico
de operarios y grúas horquilla).

Métricas oficiales del proyecto (enunciado: "eficiencia y confiabilidad", máx. 2):
  1. Confiabilidad : BER (post-FEC cuando aplica), resumida en Cobertura espacial
                     (% del área con BER ≤ 10^-6) y Disponibilidad temporal (≥ 99 %).
  2. Eficiencia    : Goodput efectivo libre de errores [Mb/s] y Eficiencia
                     Energética Incremental de Comunicación [Mbit/J].
  * Variable física intermedia: SNR [dB] para diagnóstico espacial del enlace.

Uso:
    python run_all.py          (≈ 1–2 min; resultados reproducibles con semilla 2023)
"""
from dataclasses import replace
import csv
import os
import time

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import numpy as np

from vlc_core import (
    C_LUZ,
    EMI_SOURCES,
    SNR_TH_OOK_DB,
    SNR_TH_OOK_FEC_DB,
    STATIC_OBSTACLES,
    IrUplinkParams,
    Params,
    WifiParams,
    adr_channels,
    heavy_static_boxes,
    illuminance,
    inside_any,
    ir_uplink_budget,
    los_channel,
    optimized_led_grid,
    performance,
    performance_adr,
    random_dynamic_boxes,
    receiver_grid,
    reflection_histogram,
    rms_delay_spread,
    static_boxes,
    wifi_evaluate,
    wifi_spatial_emi_dB,
)

OUT = "resultados"
os.makedirs(OUT, exist_ok=True)
plt.rcParams.update({
    "font.size": 9.5,
    "axes.grid": True,
    "grid.alpha": 0.28,
    "figure.dpi": 110,
    "savefig.dpi": 160,
    "savefig.bbox": "tight",
})

RNG = np.random.default_rng(2023)
BER_TH = 1e-6
t_start = time.time()
log_lines = []


def log(msg: str, timed: bool = True):
    line = f"[{time.time() - t_start:6.1f} s] {msg}" if timed else msg
    print(line, flush=True)
    log_lines.append(line)


def beta_of(V: float) -> float:
    return 0.0 if np.isinf(V) else 3.912 / V


def coverage(ber) -> float:
    return float(100.0 * np.mean(np.asarray(ber) <= BER_TH))


def draw_layout(ax, p: Params, boxes_dyn=None, show_leds=True, show_ap=False,
                show_static=True, show_emi=False):
    if show_static:
        for _, b in STATIC_OBSTACLES.items():
            ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2],
                                   fc="0.55", ec="k", lw=0.8, alpha=0.85))
    if boxes_dyn is not None and len(boxes_dyn):
        for b in boxes_dyn:
            c = "tab:orange" if (b[1] - b[0]) < 0.7 else "tab:red"
            ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2],
                                   fc=c, ec="k", lw=0.5))
    if show_leds:
        L = p.leds
        ax.plot(L[:, 0], L[:, 1], "*", ms=9.5, mfc="gold", mec="k", label="Luminaria VLC")
    if show_ap:
        w = WifiParams()
        ax.plot(w.ap_pos[0], w.ap_pos[1], "^", ms=8.5, mfc="tab:cyan", mec="k", label="AP Wi-Fi")
    if show_emi:
        xs_e = [v[0] for v in EMI_SOURCES.values()]
        ys_e = [v[1] for v in EMI_SOURCES.values()]
        ax.plot(xs_e, ys_e, "o", ms=7.5, mfc="none", mec="crimson", mew=1.8, label="Fuente EMI")
    ax.set_xlim(0, p.L)
    ax.set_ylim(0, p.W)
    ax.set_aspect("equal")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(False)


def write_csv(filename: str, header, rows):
    with open(f"{OUT}/{filename}", "w", newline="", encoding="utf-8-sig") as f:
        wr = csv.writer(f)
        wr.writerow(header)
        wr.writerows(rows)


csv_rows = []

# ===========================================================================
# Geometría base (nave de manufactura pesada 15 × 15 × 6 m)
# ===========================================================================
P0 = Params()
P_OPT = optimized_led_grid(P0)
SB = static_boxes()
HB = heavy_static_boxes(min_height=2.0)

rx_map, xs, ys = receiver_grid(P0, 0.25)
free_map = ~inside_any(rx_map[:, :2], SB)
# Área de tránsito activo AGV (excluye máquinas fijas y brecha muerta de 0,8 m tras estanterías de 5 m)
work_map = free_map & (rx_map[:, 1] <= 13.2)
nx, ny = len(xs), len(ys)
log(f"Grilla de mapas: {len(rx_map)} puntos ({free_map.sum()} área libre, {work_map.sum()} zona tránsito AGV)")

# ===========================================================================
# E1  Planta industrial (con fuentes EMI) e iluminancia EN 12464-1 (3×3 vs 4×4)
# ===========================================================================
E_base_empty = illuminance(P0, rx_map, None, include_reflections=True)
E_opt_empty = illuminance(P_OPT, rx_map, None, include_reflections=True)
E_base_obs = illuminance(P0, rx_map, SB, include_reflections=True)
E_opt_obs = illuminance(P_OPT, rx_map, SB, include_reflections=True)

u0_base = E_base_empty.min() / E_base_empty.mean()
u0_opt = E_opt_empty.min() / E_opt_empty.mean()
u0_base_work = E_base_obs[work_map].min() / E_base_obs[work_map].mean()
u0_opt_work = E_opt_obs[work_map].min() / E_opt_obs[work_map].mean()

fig, axs = plt.subplots(1, 3, figsize=(16.0, 4.9))
draw_layout(axs[0], P0, show_ap=True, show_emi=True)
for name, b in STATIC_OBSTACLES.items():
    axs[0].text((b[0] + b[1]) / 2, b[3] - 0.12, f"{name}\n({b[4]:.1f} m)",
                ha="center", va="top", fontsize=6.0, fontweight="bold")
axs[0].legend(loc="upper center", bbox_to_anchor=(0.5, -0.14), ncol=3, fontsize=7.8)
axs[0].set_title("Planta industrial (15 × 15 × 6 m)\nMáquinas, luminarias 3×3 y fuentes EMI")

for ax, E_arr, p_cfg, title in (
    (axs[1], E_base_empty, P0,
     f"Grilla base 3×3 (EN 12464-1: U0 < 0,60)\nmedia {E_base_empty.mean():.0f} lx · U0 = {u0_base:.2f}"),
    (axs[2], E_opt_empty, P_OPT,
     f"Grilla optimizada 4×4 (misma potencia, cumple EN)\nmedia {E_opt_empty.mean():.0f} lx · U0 = {u0_opt:.2f}"),
):
    im = ax.imshow(E_arr.reshape(ny, nx), origin="lower", extent=[0, P0.L, 0, P0.W],
                   cmap="inferno", vmin=120, vmax=450)
    cs = ax.contour(xs, ys, E_arr.reshape(ny, nx), levels=[200, 300, 400], colors="w", linewidths=0.8)
    ax.clabel(cs, fmt="%d lx", fontsize=7)
    ax.plot(p_cfg.leds[:, 0], p_cfg.leds[:, 1], "*", ms=8, mfc="gold", mec="k")
    plt.colorbar(im, ax=ax, label="Iluminancia [lx]")
    ax.set_title(title, fontsize=9.2)
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(False)

fig.savefig(f"{OUT}/E1_planta_iluminancia.png")
plt.close(fig)

csv_rows += [
    ["E1", "Grilla base 3×3 · Iluminancia media nave vacía (c/refl.) [lx]", f"{E_base_empty.mean():.0f}"],
    ["E1", "Grilla base 3×3 · Uniformidad U0 nave vacía (EN ≥ 0,60)", f"{u0_base:.2f}"],
    ["E1", "Grilla optimizada 4×4 · Iluminancia media nave vacía [lx]", f"{E_opt_empty.mean():.0f}"],
    ["E1", "Grilla optimizada 4×4 · Uniformidad U0 nave vacía (EN ≥ 0,60)", f"{u0_opt:.2f}"],
    ["E1", "Uniformidad en zona de tránsito AGV con máquinas (base / 4×4)", f"{u0_base_work:.2f} / {u0_opt_work:.2f}"],
]
log("E1 listo")

# ===========================================================================
# Canales base (mapas 0,25 m): vacío, obstáculos fijos, y fijos+móviles para S4/S5
# ===========================================================================
V_DUST = 20.0
betas_map = (0.0, beta_of(V_DUST))
h0, d0 = los_channel(P0, rx_map)
hist0, tb = reflection_histogram(P0, rx_map, None, betas_map)
log("Canal sin obstáculos listo")

dyn_map = random_dynamic_boxes(RNG, 12, 3, P0)
boxes_S45 = np.vstack([SB, dyn_map])
free_S45 = ~inside_any(rx_map[:, :2], boxes_S45)
hD, dD = los_channel(P0, rx_map, boxes_S45)
# Corrección física: en S4 y S5 las reflexiones TAMBIÉN se bloquean por obstáculos móviles
histD, _ = reflection_histogram(P0, rx_map, boxes_S45, betas_map)
log("Canal con obstáculos fijos + móviles (LOS y reflexiones) listo")

# ===========================================================================
# E2  Caso de referencia a tasa nominal industrial (10 Mb/s)
# ===========================================================================
r0 = performance(P0, h0, d0, hist0[0.0], tb)
tau = rms_delay_spread(h0, d0, hist0[0.0], tb, beta_ext=0.0)

fig, axs = plt.subplots(1, 3, figsize=(15.5, 4.4))
for ax, data, lab, cm in (
    (axs[0], 10.0 * np.log10(np.maximum(r0["P_dc"], 1e-30) * 1e3), "Potencia óptica recibida [dBm]", "viridis"),
    (axs[1], r0["snr_db"], "SNR intermedia [dB]", "viridis"),
    (axs[2], tau * 1e9, "Dispersión de retardo RMS [ns]", "magma"),
):
    im = ax.imshow(data.reshape(ny, nx), origin="lower", extent=[0, P0.L, 0, P0.W], cmap=cm)
    plt.colorbar(im, ax=ax, label=lab)
    ax.plot(P0.leds[:, 0], P0.leds[:, 1], "*", ms=8.5, mfc="gold", mec="k")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    ax.grid(False)

axs[0].set_title(f"Potencia recibida ({10*np.log10(r0['P_dc'].min()*1e3):.1f} a {10*np.log10(r0['P_dc'].max()*1e3):.1f} dBm)")
axs[1].set_title(f"SNR a 10 Mb/s ({r0['snr_db'].min():.1f} a {r0['snr_db'].max():.1f} dB)")
axs[2].set_title(f"τ_RMS (máx. {np.nanmax(tau)*1e9:.1f} ns)")
fig.suptitle("E2 · Caso de referencia (10 Mb/s): nave sin obstáculos, aire limpio, T_amb = 25 °C (Tj = 55 °C)", y=1.02)
fig.savefig(f"{OUT}/E2_referencia_potencia_snr_retardo.png")
plt.close(fig)

csv_rows += [
    ["E2", "Potencia óptica recibida min/max [dBm]", f"{10*np.log10(r0['P_dc'].min()*1e3):.1f} / {10*np.log10(r0['P_dc'].max()*1e3):.1f}"],
    ["E2", "SNR intermedia min/mediana/max a 10 Mb/s [dB]", f"{r0['snr_db'].min():.1f} / {np.median(r0['snr_db']):.1f} / {r0['snr_db'].max():.1f}"],
    ["E2", "Goodput medio caso referencia [Mb/s]", f"{r0['goodput_mbps'].mean():.2f}"],
    ["E2", "tau_RMS máximo [ns]", f"{np.nanmax(tau)*1e9:.1f}"],
]

fig, ax = plt.subplots(figsize=(7.5, 3.8))
for (x, y, c) in ((7.5, 7.5, "tab:blue"), (0.625, 0.625, "tab:red")):
    i = int(np.argmin((rx_map[:, 0] - x) ** 2 + (rx_map[:, 1] - y) ** 2))
    hr = hist0[0.0][i] / (tb[1] - tb[0])
    ax.plot(tb * 1e9, hr, c, lw=1.1, label=f"Reflexiones en ({rx_map[i,0]:.1f}, {rx_map[i,1]:.1f}) m")
    for j in range(len(P0.leds)):
        if h0[i, j] > 0:
            ax.vlines(d0[i, j] / C_LUZ * 1e9, 0, hr.max() * 1.15, colors=c, lw=1.5, alpha=0.65)
ax.set_xlim(10, 90)
ax.set_xlabel("Tiempo [ns]")
ax.set_ylabel("h(t) [1/s]")
ax.set_title("E2 · Respuesta impulsiva del canal (líneas verticales = llegadas LOS por luminaria)")
ax.legend(fontsize=8.5)
fig.savefig(f"{OUT}/E2_respuesta_impulsiva.png")
plt.close(fig)
log("E2 listo")

# ===========================================================================
# E3  Seis escenarios hostiles a 10 Mb/s (aislando cada factor y combinados)
# ===========================================================================
scen = {
    "S0 · Ideal (25 °C, limpio)":                 (P0, h0, d0, hist0[0.0], np.ones(len(rx_map), bool)),
    "S1 · Polvo V=20 m (c/scatt.) + 30 d":        (replace(P0, visibility=V_DUST, days_dirty=30), h0, d0, hist0[beta_of(V_DUST)], np.ones(len(rx_map), bool)),
    "S2 · Calor T_amb=70 °C (Tj=100 °C)":         (replace(P0, T_amb=70.0), h0, d0, hist0[0.0], np.ones(len(rx_map), bool)),
    "S3 · Luz ambiente intensa (5100 µA)":        (replace(P0, I_bg=5100e-6), h0, d0, hist0[0.0], np.ones(len(rx_map), bool)),
    "S4 · Obstáculos (máq. + 12 op. + 3 grúas)":  (P0, hD, dD, histD[0.0], free_S45),
    "S5 · Combinado (todos los factores)":        (replace(P0, visibility=V_DUST, days_dirty=30, T_amb=70.0, I_bg=5100e-6),
                                                   hD, dD, histD[beta_of(V_DUST)], free_S45),
}

res_scen = {}
fig, axs = plt.subplots(2, 3, figsize=(15.8, 9.8))
cmap = plt.get_cmap("RdYlGn")
for ax, (name, (p, h, d, hi, area)) in zip(axs.ravel(), scen.items()):
    r = performance(p, h, d, hi, tb)
    with_obs = ("Obst" in name) or ("Comb" in name)
    snr = np.where(area, r["snr_db"], np.nan)
    cov = coverage(r["ber"][area])
    gp_mean = float(np.mean(r["goodput_mbps"][area]))
    res_scen[name] = (r, cov, gp_mean, area)

    im = ax.imshow(np.clip(snr, -10, 45).reshape(ny, nx), origin="lower", extent=[0, P0.L, 0, P0.W],
                   cmap=cmap, vmin=-10, vmax=45)
    ax.contour(xs, ys, np.nan_to_num(snr, nan=99).reshape(ny, nx), levels=[SNR_TH_OOK_DB], colors="k", linewidths=1.1)
    draw_layout(ax, P0, dyn_map if with_obs else None, show_static=with_obs)
    ax.set_title(f"{name}\ncob. {cov:.1f} % · SNR med. {np.nanmedian(snr):.1f} dB · Gp {gp_mean:.2f} Mb/s", fontsize=8.8)
    csv_rows.append(["E3", f"{name} · Cobertura [%] / Goodput [Mb/s]", f"{cov:.1f} % / {gp_mean:.2f}"])
    csv_rows.append(["E3", f"{name} · SNR mediana [dB]", f"{np.nanmedian(snr):.1f}"])

fig.colorbar(im, ax=axs, shrink=0.8, label="SNR intermedia [dB]  (línea negra = 13,5 dB → BER = 10⁻⁶)")
fig.suptitle("E3 · Mapas de SNR a 10 Mb/s en seis escenarios · Cobertura = % área libre con BER ≤ 10⁻⁶\n"
             "(S4–S5: gris = máquinas, naranjo = operarios, rojo = grúas; reflexiones bloqueadas también por móviles)", y=0.99)
fig.savefig(f"{OUT}/E3_mapas_snr_escenarios.png")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.0, 4.4))
for name, (r, cov, gp_m, area) in res_scen.items():
    b = np.sort(np.log10(np.maximum(r["ber"][area], 1e-30)))
    ax.plot(b, np.arange(1, len(b) + 1) / len(b) * 100.0, lw=1.7, label=f"{name} ({cov:.1f} %)")
ax.axvline(-6.0, color="k", ls="--", lw=1.1)
ax.annotate("Umbral BER = 10⁻⁶", xy=(-6.0, 42), xytext=(-14.5, 42),
            arrowprops=dict(arrowstyle="->", color="k", lw=0.9), fontsize=8.5, fontweight="bold")
ax.set_xlim(-30, 0)
ax.set_ylim(0, 103)
ax.set_xlabel("log10(BER)")
ax.set_ylabel("% del área libre con BER ≤ x")
ax.set_title("E3 · Distribución acumulada (CDF) de BER por escenario (10 Mb/s)")
ax.legend(fontsize=7.8, loc="lower left", framealpha=0.95)
fig.savefig(f"{OUT}/E3_cdf_ber.png")
plt.close(fig)
log("E3 listo")

# ===========================================================================
# Grilla gruesa para barridos paramétricos (0,5 m) en zona libre
# ===========================================================================
rx_s, _, _ = receiver_grid(P0, 0.5)
free_s = ~inside_any(rx_s[:, :2], SB)
work_s = free_s & (rx_s[:, 1] <= 13.2)
rx_sf = rx_s[free_s]
rx_sw = rx_s[work_s]

V_list = [np.inf, 200.0, 100.0, 50.0, 30.0, 20.0, 15.0, 10.0, 7.0, 5.0]
betas_s = tuple(sorted(set(beta_of(v) for v in V_list)))
h_s, d_s = los_channel(P0, rx_sf, SB)
hist_s, tb = reflection_histogram(P0, rx_sf, SB, betas_s)
h_s0, d_s0 = los_channel(P0, rx_sf)
hist_s0, _ = reflection_histogram(P0, rx_sf, None, betas_s)
log(f"Grilla de barridos lista ({len(rx_sf)} pts libres, {len(rx_sw)} pts tránsito AGV)")

# ===========================================================================
# E4  Polvo en suspensión (Beer-Lambert + scattering Henyey-Greenstein)
# ===========================================================================
fig, axs = plt.subplots(1, 2, figsize=(12.4, 4.4))
Vplot = [v for v in V_list if not np.isinf(v)]
rows_e4 = []
for Rb, blue, lab_r, c in (
    (2e6, False, "2 Mb/s (base)", "tab:green"),
    (10e6, False, "10 Mb/s (nominal)", "tab:blue"),
    (20e6, True, "20 Mb/s (c/filtro azul)", "tab:red"),
):
    med, p95, cov_l = [], [], []
    for v in Vplot:
        p = replace(P0, visibility=v, Rb=Rb, use_blue_filter=blue)
        r = performance(p, h_s0, d_s0, hist_s0[beta_of(v)], tb)
        lb = np.log10(np.maximum(r["ber"], 1e-40))
        med.append(np.median(lb))
        p95.append(np.percentile(lb, 95))
        cov_l.append(coverage(r["ber"]))
        rows_e4.append([
            Rb / 1e6, v, f"{np.median(r['snr_db']):.2f}", f"{np.percentile(r['snr_db'], 5):.2f}",
            f"{coverage(r['ber']):.1f}", f"{np.mean(r['goodput_mbps']):.2f}",
        ])
    axs[0].semilogx(Vplot, med, "o-", c=c, label=f"{lab_r} · mediana")
    axs[0].semilogx(Vplot, p95, "s--", c=c, alpha=0.65, label=f"{lab_r} · peor 5 %")
    axs[1].semilogx(Vplot, cov_l, "o-", c=c, label=lab_r)

axs[0].axhline(-6, color="k", ls=":", lw=1.0)
axs[0].set_ylim(-30, 0)
axs[0].set_xlabel("Visibilidad V [m] (menor = más polvo)")
axs[0].set_ylabel("log10(BER)")
axs[0].set_title("E4 · BER vs polvo en suspensión (extinción + scattering)")
axs[0].legend(fontsize=7.2, ncol=2, loc="lower left")

axs[1].set_xlabel("Visibilidad V [m] (menor = más polvo)")
axs[1].set_ylabel("Cobertura con BER ≤ 10⁻⁶ [%]")
axs[1].set_ylim(-2, 105)
axs[1].set_title("E4 · Cobertura confiable vs concentración de polvo")
axs[1].legend(fontsize=8.2, loc="lower left")

for ax in axs:
    ax.invert_xaxis()
    ax.set_xticks([200, 100, 50, 20, 10, 5])
    ax.set_xticklabels(["200", "100", "50", "20", "10", "5"])

fig.savefig(f"{OUT}/E4_polvo_ber.png")
plt.close(fig)
write_csv("E4_polvo.csv", ["Rb_Mbps", "visibilidad_m", "SNR_mediana_dB", "SNR_p5_dB", "cobertura_%", "goodput_medio_Mbps"], rows_e4)
log("E4 listo")

# ===========================================================================
# E5  Suciedad en ópticas (LDD + acoplamiento térmico Tj y luz ambiente I_bg)
# ===========================================================================
days = np.arange(0, 121, 5)
fig, ax = plt.subplots(figsize=(7.8, 4.3))
rows_e5 = []
for v, c in ((np.inf, "tab:green"), (50.0, "tab:blue"), (20.0, "tab:orange"), (10.0, "tab:red")):
    cov_l = []
    for dd in days:
        p = replace(P0, visibility=v, days_dirty=float(dd))
        r = performance(p, h_s, d_s, hist_s[beta_of(v)], tb)
        cv = coverage(r["ber"])
        cov_l.append(cv)
        rows_e5.append([v, int(dd), f"{cv:.1f}", f"{np.median(r['snr_db']):.2f}", f"{np.mean(r['goodput_mbps']):.2f}"])
    lab = "Aire limpio" if np.isinf(v) else f"Polvo V = {v:.0f} m (c/scattering)"
    ax.plot(days, cov_l, "o-", ms=3.2, c=c, label=lab)
ax.set_xlabel("Días sin limpiar ópticas de luminarias y receptores AGV")
ax.set_ylabel("Cobertura con BER ≤ 10⁻⁶ [%]")
ax.set_ylim(-2, 102)
ax.set_title("E5 · Efecto de suciedad acumulada LDD + acople térmico/I_bg (con máquinas, 10 Mb/s)")
ax.legend(fontsize=8.2)
fig.savefig(f"{OUT}/E5_suciedad_limpieza.png")
plt.close(fig)
write_csv("E5_suciedad.csv", ["visibilidad_m", "dias", "cobertura_%", "SNR_mediana_dB", "goodput_medio_Mbps"], rows_e5)
log("E5 listo")

# ===========================================================================
# E6  Temperatura ambiente y temperatura de unión Tj del LED (Davis et al. 2019)
# ===========================================================================
Ts = np.arange(0, 91, 5)
fig, axs = plt.subplots(1, 2, figsize=(12.6, 4.5))

# Régimen realista limitado por ruido (polvo V=15 m + 30 días suciedad + luz solar/horno)
# donde la caída térmica del flujo del LED impacta directamente la SNR (en aire limpio la ISI domina)
p_th_ref = replace(P0, visibility=15.0, days_dirty=30.0, I_bg=5100e-6)
ref_t = performance(p_th_ref, h_s0, d_s0, hist_s0[beta_of(15.0)], tb)
ref_snr = float(np.median(ref_t["snr_db"]))
ref_clean = float(np.median(performance(P0, h_s0, d_s0, hist_s0[0.0], tb)["snr_db"]))

for kT, c, ls, lab in (
    (0.002, "tab:green", "-", "Droop 0,2 %/°C (sin derating driver)"),
    (0.004, "tab:blue", "-", "Droop 0,4 %/°C (sin derating driver)"),
    (0.010, "tab:red", "-", "Droop 1,0 %/°C (sin derating driver)"),
):
    dsnr = [
        float(np.median(performance(replace(p_th_ref, T_amb=float(T), k_T_led=kT, T_j_derate=None),
                                    h_s0, d_s0, hist_s0[beta_of(15.0)], tb)["snr_db"])) - ref_snr
        for T in Ts
    ]
    axs[0].plot(Ts, dsnr, ls, c=c, lw=1.6, label=lab)

dsnr_prot = [
    float(np.median(performance(replace(p_th_ref, T_amb=float(T), k_T_led=0.004, T_j_derate=85.0),
                                h_s0, d_s0, hist_s0[beta_of(15.0)], tb)["snr_db"])) - ref_snr
    for T in Ts
]
dsnr_clean = [
    float(np.median(performance(replace(P0, T_amb=float(T), k_T_led=0.004, T_j_derate=85.0),
                                h_s0, d_s0, hist_s0[0.0], tb)["snr_db"])) - ref_clean
    for T in Ts
]
axs[0].plot(Ts, dsnr_prot, "--", c="k", lw=2.0, label="0,4 %/°C + derating driver (Tj > 85 °C, hostil)")
axs[0].plot(Ts, dsnr_clean, ":", c="tab:purple", lw=1.8, label="Ídem en aire limpio (régimen limitado por ISI)")
axs[0].axvline(53.0, color="0.45", ls=":", lw=1.0)
axs[0].annotate("Tj ≈ 85 °C\n(inicio derating)", xy=(53.0, dsnr_prot[11]), xytext=(26.0, min(dsnr_prot) * 0.65),
                arrowprops=dict(arrowstyle="->", color="0.3", lw=0.9), fontsize=7.8, color="0.25")
axs[0].set_xlabel("Temperatura ambiente T_amb [°C]  (Tj ≈ T_amb + 30–32 °C)")
axs[0].set_ylabel("Variación de SNR mediana respecto a 25 °C [dB]")
axs[0].set_title("E6 · Pérdida térmica de SNR (V=15 m, 30 d, 5100 µA vs aire limpio)")
axs[0].legend(fontsize=7.3, loc="lower left")

# Descomposición de varianzas de ruido en el centro de la nave (10 Mb/s)
i_ctr = int(np.argmin((rx_sf[:, 0] - 7.25) ** 2 + (rx_sf[:, 1] - 6.25) ** 2))
labels_t, sh_l, th_l, dk_l = [], [], [], []
for T in (25, 50, 70, 85):
    p_t = replace(P0, T_amb=float(T))
    r_t = performance(p_t, h_s0[i_ctr:i_ctr + 1], d_s0[i_ctr:i_ctr + 1], hist_s0[0.0][i_ctr:i_ctr + 1], tb)
    I_dark = p_t.I_dark0 * (2.0 ** ((T - 25.0) / 10.0))
    s_dark = 2.0 * 1.602176634e-19 * I_dark * p_t.Rb
    labels_t.append(f"{T} °C\n(Tj={p_t.T_j:.0f}°C)")
    sh_l.append(max(float(r_t["s_shot"][0]) - s_dark, 1e-25))
    th_l.append(float(r_t["s_th"]))
    dk_l.append(s_dark)

x_t = np.arange(len(labels_t))
axs[1].bar(x_t - 0.27, sh_l, 0.27, label="Shot (señal + luz ambiente)")
axs[1].bar(x_t, th_l, 0.27, label="Térmico preamplificador")
axs[1].bar(x_t + 0.27, dk_l, 0.27, label="Corriente oscura")
axs[1].set_yscale("log")
axs[1].set_xticks(x_t)
axs[1].set_xticklabels(labels_t)
axs[1].set_ylabel("Varianza de ruido [A²]")
axs[1].set_title("E6 · Componentes de ruido vs temperatura (centro, 10 Mb/s)")
axs[1].legend(fontsize=7.8, loc="upper center", bbox_to_anchor=(0.5, -0.15), ncol=3)

fig.savefig(f"{OUT}/E6_temperatura.png")
plt.close(fig)

idx70 = int(np.where(Ts == 70)[0][0])
idx85 = int(np.where(Ts == 85)[0][0])
csv_rows.append(["E6", "ΔSNR mediana a T_amb=70 °C (Tj=102 °C, régimen ruido hostil) [dB]", f"{dsnr_prot[idx70]:.2f}"])
csv_rows.append(["E6", "ΔSNR mediana a T_amb=85 °C (Tj=117 °C, régimen ruido hostil) [dB]", f"{dsnr_prot[idx85]:.2f}"])
log("E6 listo")

# ===========================================================================
# E7  Tasa de bits, ancho de banda del LED (Le Minh 2009), ISI y Goodput
# ===========================================================================
Rbs = np.array([1, 2, 5, 10, 15, 20, 30, 50, 75, 100]) * 1e6
fig, axs = plt.subplots(1, 2, figsize=(12.8, 4.5))
rows_e7 = []
cases_e7 = (
    ("LED crudo 3 MHz (Ideal)",            replace(P0, f_3db_led=3e6, use_blue_filter=False), h_s0, d_s0, hist_s0[0.0], "tab:gray", "--"),
    ("Base pre-ec. 15 MHz (Ideal)",        replace(P0, f_3db_led=15e6, use_blue_filter=False), h_s0, d_s0, hist_s0[0.0], "tab:green", "-"),
    ("Filtro azul 25 MHz (Ideal)",         replace(P0, use_blue_filter=True), h_s0, d_s0, hist_s0[0.0], "tab:blue", "-"),
    ("Combinado (Base 15 MHz, sin FEC)",   replace(P0, visibility=20.0, days_dirty=30.0, T_amb=70.0, I_bg=5100e-6), h_s, d_s, hist_s[beta_of(20.0)], "tab:red", "-"),
    ("Combinado (Filtro azul 25 MHz+FEC)", replace(P0, visibility=20.0, days_dirty=7.0, T_amb=70.0, I_bg=5100e-6, use_blue_filter=True, use_fec=True), h_s, d_s, hist_s[beta_of(20.0)], "tab:purple", "-"),
)
for name, p_cfg, h, d, hi, c, ls in cases_e7:
    cov_l, gp_l = [], []
    for Rb in Rbs:
        r = performance(replace(p_cfg, Rb=float(Rb)), h, d, hi, tb)
        cv = coverage(r["ber"])
        gm = float(np.mean(r["goodput_mbps"]))
        cov_l.append(cv)
        gp_l.append(gm)
        rows_e7.append([name, Rb / 1e6, f"{cv:.1f}", f"{np.median(r['snr_db']):.2f}", f"{gm:.2f}"])
    axs[0].semilogx(Rbs / 1e6, cov_l, "o" + ls, ms=4, c=c, label=name)
    axs[1].semilogx(Rbs / 1e6, gp_l, "o" + ls, ms=4, c=c, label=name)

axs[0].set_xlabel("Tasa de bits OOK [Mb/s]")
axs[0].set_ylabel("Cobertura con BER ≤ 10⁻⁶ [%]")
axs[0].set_ylim(-2, 104)
axs[0].set_title("E7a · Confiabilidad vs tasa (efecto del BW del LED e ISI)")
axs[0].legend(fontsize=7.4, loc="lower left")

axs[1].set_xlabel("Tasa de bits OOK [Mb/s]")
axs[1].set_ylabel("Goodput efectivo medio [Mb/s]")
axs[1].set_title("E7b · Eficiencia (Goodput útil) vs tasa de transmisión")
axs[1].legend(fontsize=7.4, loc="upper left")

fig.savefig(f"{OUT}/E7_tasa_bits.png")
plt.close(fig)
write_csv("E7_tasa_bits.csv", ["escenario", "Rb_Mbps", "cobertura_%", "SNR_mediana_dB", "goodput_medio_Mbps"], rows_e7)
log("E7 listo")

# ===========================================================================
# E8  Obstáculos móviles (Monte Carlo equiparado a 10 Mb/s): 1-PD, ADR, Wi-Fi e Híbrido
# ===========================================================================
N_MC = 150
rx_mc = rx_sw[RNG.choice(len(rx_sw), 180, replace=False)]
W_PLANT = WifiParams(Rb_target=10e6, use_fec=False)
SF_MC = RNG.standard_normal(len(rx_mc))

# Precalcular reflexiones estáticas para cada configuración VLC
p_9_60 = replace(P0, fov_deg=60.0)
p_9_80 = replace(P0, fov_deg=80.0)
p_4_60 = replace(P0, led_xy=((3.75, 3.75), (3.75, 11.25), (11.25, 3.75), (11.25, 11.25)),
                 P_opt=45.0 * 9.0 / 4.0)
p_9_adr = replace(P0, adr_enabled=True)
p_16_adr = replace(P_OPT, adr_enabled=True)

hist_9_60, _ = reflection_histogram(p_9_60, rx_mc, SB, (0.0,))
hist_9_80, _ = reflection_histogram(p_9_80, rx_mc, SB, (0.0,))
hist_4_60, _ = reflection_histogram(p_4_60, rx_mc, SB, (0.0,))
_, _, hist_9_adr_l, _ = adr_channels(p_9_adr, rx_mc, SB, (0.0,), include_refl=True)
hist_9_adr = [hb[0.0] for hb in hist_9_adr_l]
_, _, hist_16_adr_l, _ = adr_channels(p_16_adr, rx_mc, SB, (0.0,), include_refl=True)
hist_16_adr = [hb[0.0] for hb in hist_16_adr_l]

n_workers_list = [0, 5, 10, 20, 30, 40]
e8_names = [
    "VLC · 9 lum., 1-PD FOV 60°",
    "VLC · 9 lum., 1-PD FOV 80°",
    "VLC · 4 lum., 1-PD FOV 60°",
    "VLC · 9 lum., ADR piramidal",
    "VLC · 16 lum. (EN12464) + ADR",
    "Wi-Fi 2,4 GHz (EMI espacial)",
    "Híbrido VLC (9 lum. ADR) + Wi-Fi",
]
avail_e8 = {k: [] for k in e8_names}
gp_e8 = {k: [] for k in e8_names}

for nw in n_workers_list:
    ok = {k: np.zeros(len(rx_mc)) for k in e8_names}
    val_cnt = np.zeros(len(rx_mc))
    gp_acc = {k: [] for k in e8_names}
    for _ in range(N_MC):
        dyn = random_dynamic_boxes(RNG, nw, nw // 10, P0)
        # Excluir solapamiento físico ilegal (un operario/grúa no puede ocupar la misma coordenada exacta del AGV)
        valid = ~inside_any(rx_mc[:, :2], dyn)
        val_cnt += valid
        boxes = np.vstack([SB, dyn]) if len(dyn) else SB

        # 1-PD configuraciones
        for nm, p_c, hi in (
            ("VLC · 9 lum., 1-PD FOV 60°", p_9_60, hist_9_60[0.0]),
            ("VLC · 9 lum., 1-PD FOV 80°", p_9_80, hist_9_80[0.0]),
            ("VLC · 4 lum., 1-PD FOV 60°", p_4_60, hist_4_60[0.0]),
        ):
            h_i, d_i = los_channel(p_c, rx_mc, boxes)
            r_i = performance(p_c, h_i, d_i, hi, tb)
            ok[nm] += (r_i["ber"] <= BER_TH) & valid
            gp_acc[nm].append(float(np.mean(r_i["goodput_mbps"][valid])))

        # ADR 9 y 16 luminarias
        ha9, da9, _, _ = adr_channels(p_9_adr, rx_mc, boxes, include_refl=False)
        r_a9 = performance_adr(p_9_adr, ha9, da9, hist_9_adr, tb)
        ok["VLC · 9 lum., ADR piramidal"] += (r_a9["ber"] <= BER_TH) & valid
        gp_acc["VLC · 9 lum., ADR piramidal"].append(float(np.mean(r_a9["goodput_mbps"][valid])))

        ha16, da16, _, _ = adr_channels(p_16_adr, rx_mc, boxes, include_refl=False)
        r_a16 = performance_adr(p_16_adr, ha16, da16, hist_16_adr, tb)
        ok["VLC · 16 lum. (EN12464) + ADR"] += (r_a16["ber"] <= BER_TH) & valid
        gp_acc["VLC · 16 lum. (EN12464) + ADR"].append(float(np.mean(r_a16["goodput_mbps"][valid])))

        # Wi-Fi 2.4 GHz (misma tasa 10 Mb/s, EMI espacial Middleton Class A)
        rw = wifi_evaluate(W_PLANT, rx_mc, HB, dyn, SF_MC, RNG, n_fad=1, use_spatial_emi=True)
        w_ok = rw["ber"] <= BER_TH
        ok["Wi-Fi 2,4 GHz (EMI espacial)"] += w_ok & valid
        gp_acc["Wi-Fi 2,4 GHz (EMI espacial)"].append(float(np.mean(rw["goodput_mbps"][valid])))

        # Híbrido
        ok["Híbrido VLC (9 lum. ADR) + Wi-Fi"] += ((r_a9["ber"] <= BER_TH) | w_ok) & valid
        gp_acc["Híbrido VLC (9 lum. ADR) + Wi-Fi"].append(
            float(np.mean(np.maximum(r_a9["goodput_mbps"][valid], rw["goodput_mbps"][valid])))
        )

    for k in e8_names:
        av = ok[k] / np.maximum(val_cnt, 1.0)
        avail_e8[k].append(float(100.0 * np.mean(av >= 0.99)))
        gp_e8[k].append(float(np.mean(gp_acc[k])))
    log(f"E8: {nw} operarios listo")

fig, ax = plt.subplots(figsize=(9.2, 4.9))
styles_e8 = ["o-", "s--", "v:", "D-", "^-", "X--", "P-"]
colors_e8 = ["tab:blue", "tab:green", "tab:gray", "tab:orange", "tab:purple", "tab:cyan", "k"]
for k, st, c in zip(e8_names, styles_e8, colors_e8):
    ax.plot(n_workers_list, avail_e8[k], st, c=c, lw=1.7, ms=5.5, label=k)
ax.set_xlabel("Operarios en la nave (+1 grúa horquilla por cada 10)")
ax.set_ylabel("% de ubicaciones con disponibilidad ≥ 99 %")
ax.set_ylim(0, 104)
ax.set_title(f"E8 · Bloqueo por obstáculos móviles a 10 Mb/s (Monte Carlo equiparado, {N_MC} realizaciones)")
ax.legend(fontsize=7.8, loc="lower left", ncol=2, framealpha=0.93)
fig.savefig(f"{OUT}/E8_obstaculos_montecarlo.png")
plt.close(fig)

write_csv(
    "E8_obstaculos.csv",
    ["configuracion"] + [f"{n}_operarios_avail99%" for n in n_workers_list],
    [[k] + [f"{x:.1f}" for x in avail_e8[k]] for k in e8_names],
)
log("E8 listo")

# ===========================================================================
# E9  Comparación equiparada VLC vs Wi-Fi: EMI Middleton Class A, Polvo, Mapa EMI y Eficiencia
# ===========================================================================
# Mismo conjunto de 150 realizaciones Monte Carlo (10 operarios + 1 grúa) para VLC y Wi-Fi a 10 Mb/s
N_MC_E9 = 150
dyn_e9 = [random_dynamic_boxes(RNG, 10, 1, P0) for _ in range(N_MC_E9)]
val_e9 = np.array([~inside_any(rx_mc[:, :2], d) for d in dyn_e9])     # (N_MC_E9, Nrx)
val_sum_e9 = np.maximum(val_e9.sum(axis=0), 1.0)

# Precalcular LOS VLC para las 200 realizaciones (9 lum 1-PD y 16 lum ADR)
h9_e9, d9_e9 = [], []
ha16_e9, da16_e9 = [], []
for d_dyn in dyn_e9:
    bx = np.vstack([SB, d_dyn])
    h_, d_ = los_channel(p_9_60, rx_mc, bx)
    h9_e9.append(h_)
    d9_e9.append(d_)
    ha_, da_, _, _ = adr_channels(p_16_adr, rx_mc, bx, include_refl=False)
    ha16_e9.append(ha_)
    da16_e9.append(da_)

# Reflexiones para aire limpio (0,0) y cada visibilidad de Vplot en E9b
betas_e9 = tuple(sorted(set([0.0] + [beta_of(v) for v in Vplot])))
hist9_e9_v, _ = reflection_histogram(p_9_60, rx_mc, SB, betas_e9)
_, _, hist16_e9_v_l, _ = adr_channels(p_16_adr, rx_mc, SB, betas_e9, include_refl=True)


def eval_vlc_mc_e9(p_cfg, is_adr: bool, vis: float):
    b_ext = beta_of(vis)
    p_run = replace(p_cfg, visibility=vis)
    ok_c = np.zeros(len(rx_mc))
    gps = []
    for it in range(N_MC_E9):
        if is_adr:
            hi = [hb[b_ext] for hb in hist16_e9_v_l]
            r = performance_adr(p_run, ha16_e9[it], da16_e9[it], hi, tb)
        else:
            r = performance(p_run, h9_e9[it], d9_e9[it], hist9_e9_v[b_ext], tb)
        ok_c += (r["ber"] <= BER_TH) & val_e9[it]
        gps.append(float(np.mean(r["goodput_mbps"][val_e9[it]])))
    av = ok_c / val_sum_e9
    return float(100.0 * np.mean(av >= 0.99)), float(np.mean(gps))


def eval_wifi_mc_e9(w_cfg: WifiParams, emi_override=None, use_spatial=False):
    ok_c = np.zeros(len(rx_mc))
    gps = []
    snrs = []
    for it in range(N_MC_E9):
        rw = wifi_evaluate(w_cfg, rx_mc, HB, dyn_e9[it], SF_MC, RNG, n_fad=1,
                           emi_dB_override=emi_override, use_spatial_emi=use_spatial)
        ok_c += (rw["ber"] <= BER_TH) & val_e9[it]
        gps.append(float(np.mean(rw["goodput_mbps"][val_e9[it]])))
        snrs.append(float(np.median(rw["snr_db"][val_e9[it]])))
    av = ok_c / val_sum_e9
    return float(100.0 * np.mean(av >= 0.99)), float(np.mean(gps)), float(np.mean(snrs))


emis_e9 = [0, 10, 15, 20, 25, 30, 35, 40]
wifi_emi_raw = [eval_wifi_mc_e9(W_PLANT, emi_override=float(e)) for e in emis_e9]
wifi_emi_fec = [eval_wifi_mc_e9(replace(W_PLANT, use_fec=True), emi_override=float(e)) for e in emis_e9]
vlc9_clean_av, vlc9_clean_gp = eval_vlc_mc_e9(p_9_60, False, np.inf)
vlc16_clean_av, vlc16_clean_gp = eval_vlc_mc_e9(replace(p_16_adr, use_fec=True), True, np.inf)

vlc9_v_res = [eval_vlc_mc_e9(p_9_60, False, v) for v in Vplot]
vlc16_v_res = [eval_vlc_mc_e9(replace(p_16_adr, use_fec=True), True, v) for v in Vplot]
wifi_sp_raw_av, wifi_sp_raw_gp, _ = eval_wifi_mc_e9(W_PLANT, use_spatial=True)
wifi_sp_fec_av, wifi_sp_fec_gp, _ = eval_wifi_mc_e9(replace(W_PLANT, use_fec=True), use_spatial=True)

fig, axs = plt.subplots(1, 2, figsize=(12.6, 4.5))
axs[0].plot(emis_e9, [x[0] for x in wifi_emi_raw], "D--", c="tab:cyan", label="Wi-Fi 2,4 GHz (10 Mb/s, sin FEC)")
axs[0].plot(emis_e9, [x[0] for x in wifi_emi_fec], "s-", c="teal", label="Wi-Fi 2,4 GHz (10 Mb/s + RS-FEC)")
axs[0].plot(emis_e9, [vlc9_clean_av] * len(emis_e9), "o--", c="tab:blue", label="VLC base 9 lum. (10 Mb/s, inmune EMI)")
axs[0].plot(emis_e9, [vlc16_clean_av] * len(emis_e9), "^-", c="tab:purple", label="VLC 16 lum. ADR + RS-FEC (inmune EMI)")
for x_e, lbl_e in ((5, "Oficina"), (15, "Planta\nligera"), (24, "VFD /\nmotores"), (34, "Horno /\nsoldadura")):
    axs[0].text(x_e, 4, lbl_e, ha="center", fontsize=7.0, color="0.35")
axs[0].set_xlabel("Severidad media de EMI impulsiva Middleton Class A (A=0,1) [dB]")
axs[0].set_ylabel("% ubicaciones con disponibilidad ≥ 99 %")
axs[0].set_ylim(0, 105)
axs[0].set_title("E9a · Interferencia electromagnética industrial (equiparado a 10 Mb/s)")
axs[0].legend(fontsize=7.5, loc="center left")

axs[1].semilogx(Vplot, [x[0] for x in vlc9_v_res], "o--", c="tab:blue", label="VLC base 9 lum. (10 Mb/s)")
axs[1].semilogx(Vplot, [x[0] for x in vlc16_v_res], "^-", c="tab:purple", label="VLC 16 lum. ADR + RS-FEC")
axs[1].semilogx(Vplot, [wifi_sp_raw_av] * len(Vplot), "D--", c="tab:cyan", label="Wi-Fi (EMI espacial planta, sin FEC)")
axs[1].semilogx(Vplot, [wifi_sp_fec_av] * len(Vplot), "s-", c="teal", label="Wi-Fi (EMI espacial planta + RS-FEC)")
axs[1].invert_xaxis()
axs[1].set_xticks([200, 100, 50, 20, 10, 5])
axs[1].set_xticklabels(["200", "100", "50", "20", "10", "5"])
axs[1].set_xlabel("Visibilidad V [m] (menor = más polvo, con scattering)")
axs[1].set_ylabel("% ubicaciones con disponibilidad ≥ 99 %")
axs[1].set_ylim(0, 105)
axs[1].set_title("E9b · Polvo en suspensión (complementariedad VLC vs Wi-Fi)")
axs[1].legend(fontsize=7.5, loc="lower left")

fig.savefig(f"{OUT}/E9_comparacion_wifi.png")
plt.close(fig)

write_csv(
    "E9_wifi.csv",
    ["EMI_dB", "wifi_sinFEC_avail99%", "wifi_FEC_avail99%", "wifi_FEC_goodput_Mbps",
     "vlc_base_avail99%", "vlc_ADR_FEC_avail99%"],
    [
        [e, f"{r0_[0]:.1f}", f"{rf_[0]:.1f}", f"{rf_[1]:.2f}", f"{vlc9_clean_av:.1f}", f"{vlc16_clean_av:.1f}"]
        for e, r0_, rf_ in zip(emis_e9, wifi_emi_raw, wifi_emi_fec)
    ],
)

# --- E9c: Mapa espacial de EMI, mapa de disponibilidad Wi-Fi bajo EMI y eficiencia energética + Uplink IR ---
emi_map_full = wifi_spatial_emi_dB(rx_map, W_PLANT)
SF_full = RNG.standard_normal(len(rx_map))
rw_map = wifi_evaluate(replace(W_PLANT, use_fec=True), rx_map, HB, None, SF_full, RNG, n_fad=200, use_spatial_emi=True)
wifi_av_map = 100.0 * np.mean(rw_map["ber"] <= BER_TH, axis=1)
wifi_av_map_masked = np.where(free_map, wifi_av_map, np.nan)

# Verificación de Uplink IR 850 nm (bidireccional libre de RF) en el escenario hostil con grilla 4×4
ir_res_base = ir_uplink_budget(replace(P0, visibility=20.0, days_dirty=7.0, T_amb=70.0), rx_sw, SB)
ir_res_opt = ir_uplink_budget(replace(P_OPT, visibility=20.0, days_dirty=7.0, T_amb=70.0), rx_sw, SB)
ir_cov_opt = coverage(ir_res_opt["ber"])
ir_gp_opt = float(np.mean(ir_res_opt["goodput_mbps"]))

# Eficiencia energética incremental [Mbit/J] (en planta con EMI espacial, aire limpio y 10 op + 1 grúa)
p_16_adr_fec = replace(p_16_adr, use_fec=True)
ee_labels = ["VLC base\n(9 lum.)", "VLC ADR+4×4\n+RS-FEC", "Wi-Fi 2,4 GHz\n(EMI planta)", "Wi-Fi + FEC\n(EMI planta)"]
gp_vals = [vlc9_clean_gp, vlc16_clean_gp, wifi_sp_raw_gp, wifi_sp_fec_gp]
ee_vals = [
    vlc9_clean_gp / p_9_60.P_comm_total_W,
    vlc16_clean_gp / p_16_adr_fec.P_comm_total_W,
    wifi_sp_raw_gp / W_PLANT.P_comm_total_W,
    wifi_sp_fec_gp / W_PLANT.P_comm_total_W,
]

fig, axs = plt.subplots(1, 3, figsize=(16.2, 4.7))
im0 = axs[0].imshow(emi_map_full.reshape(ny, nx), origin="lower", extent=[0, P0.L, 0, P0.W], cmap="hot_r")
draw_layout(axs[0], P0, show_leds=False, show_ap=True, show_static=True, show_emi=True)
plt.colorbar(im0, ax=axs[0], label="EMI sobre ruido térmico [dB]")
axs[0].set_title(f"E9c · Mapa espacial de EMI industrial\n({emi_map_full.min():.0f} a {emi_map_full.max():.0f} dB; picos en horno/soldadura)")

im1 = axs[1].imshow(wifi_av_map_masked.reshape(ny, nx), origin="lower", extent=[0, P0.L, 0, P0.W],
                    cmap="RdYlGn", vmin=50, vmax=100)
draw_layout(axs[1], P0, show_leds=False, show_ap=True, show_static=True, show_emi=True)
plt.colorbar(im1, ax=axs[1], label="Disponibilidad Wi-Fi + FEC [%]")
axs[1].set_title(f"Disponibilidad Wi-Fi bajo EMI real de nave\n(ubicaciones ≥ 99 %: {wifi_sp_fec_av:.1f} %)")

x_ee = np.arange(len(ee_labels))
w_bar = 0.36
b1 = axs[2].bar(x_ee - w_bar / 2, gp_vals, w_bar, color="tab:blue", label="Goodput medio [Mb/s]")
ax_r = axs[2].twinx()
b2 = ax_r.bar(x_ee + w_bar / 2, ee_vals, w_bar, color="tab:green", label="Eficiencia energ. [Mbit/J]")
axs[2].set_xticks(x_ee)
axs[2].set_xticklabels(ee_labels, fontsize=7.8)
axs[2].set_ylabel("Goodput efectivo [Mb/s]", color="tab:blue")
ax_r.set_ylabel("Eficiencia incremental [Mbit/J]", color="tab:green")
axs[2].set_ylim(0, 12.0)
ax_r.set_ylim(0, 1.25)
for bar, v in zip(b1, gp_vals):
    axs[2].text(bar.get_x() + bar.get_width() / 2, v + 0.2, f"{v:.1f}", ha="center", fontsize=7.2, color="tab:blue")
for bar, v in zip(b2, ee_vals):
    ax_r.text(bar.get_x() + bar.get_width() / 2, v + 0.025, f"{v:.2f}", ha="center", fontsize=7.2, color="tab:green")
axs[2].set_title(f"Eficiencia (Goodput y Mbit/J) a 10 Mb/s\nUplink IR 850 nm: cob. {ir_cov_opt:.1f} % ({ir_gp_opt:.2f} Mb/s)")
lines1, labels1 = axs[2].get_legend_handles_labels()
lines2, labels2 = ax_r.get_legend_handles_labels()
axs[2].legend(lines1 + lines2, labels1 + labels2, fontsize=7.3, loc="upper right")
ax_r.grid(False)

fig.savefig(f"{OUT}/E9_mapa_espacial_eficiencia.png")
plt.close(fig)

csv_rows += [
    ["E9", "Wi-Fi (EMI espacial planta, sin FEC / +FEC) · ubicaciones ≥ 99 % [%]", f"{wifi_sp_raw_av:.1f} / {wifi_sp_fec_av:.1f}"],
    ["E9", "VLC (base 9 lum. / 16 lum. ADR+FEC) · ubicaciones ≥ 99 % [%]", f"{vlc9_clean_av:.1f} / {vlc16_clean_av:.1f}"],
    ["E9", "Eficiencia energética incremental VLC ADR+4×4 vs Wi-Fi+FEC [Mbit/J]", f"{ee_vals[1]:.2f} vs {ee_vals[3]:.2f}"],
    ["E9", "Uplink IR 850 nm (4×4, V=20 m, 7 d, 70 °C) · Cobertura [%] / Goodput [Mb/s]", f"{ir_cov_opt:.1f} % / {ir_gp_opt:.2f}"],
]
log("E9 listo")

# ===========================================================================
# E10  Estudio de ablación del diseño mejorado (escenario combinado hostil, Monte Carlo)
# ===========================================================================
hard = dict(visibility=20.0, T_amb=70.0, I_bg=5100e-6)
d_base = replace(P0, days_dirty=30.0, **hard)
d_filt = replace(P0, days_dirty=7.0, use_blue_filter=True, **hard)
d_fec = replace(d_filt, use_fec=True)
d_auto = replace(optimized_led_grid(d_fec), adr_enabled=True)
w_e10 = replace(W_PLANT, use_fec=True)

beta20 = d_base.beta_ext
hist_e10_base, _ = reflection_histogram(d_base, rx_mc, SB, (beta20,))
_, _, hist_e10_adr_l, _ = adr_channels(d_auto, rx_mc, SB, (beta20,), include_refl=True)
hist_e10_adr = [hb[beta20] for hb in hist_e10_adr_l]

e10_keys = [
    "1. Base\nVLC",
    "2. +Filtro\n+7 días",
    "3. +RS\nFEC",
    "4. +ADR+4×4\nAutónomo",
    "5. Wi-Fi\n+FEC",
    "6. Híbrido\nVLC+Wi-Fi",
]
okD = {k: np.zeros(len(rx_mc)) for k in e10_keys}
gpD = {k: [] for k in e10_keys}
eeD = {k: [] for k in e10_keys}
val_e10 = np.zeros(len(rx_mc))
snr_cdf_keys = ("1. Base\nVLC", "2. +Filtro\n+7 días", "4. +ADR+4×4\nAutónomo")
snr_store = {k: [] for k in snr_cdf_keys}

for _ in range(N_MC):
    dyn = random_dynamic_boxes(RNG, 12, 3, P0)
    valid = ~inside_any(rx_mc[:, :2], dyn)
    val_e10 += valid
    boxes = np.vstack([SB, dyn])

    h1, d1 = los_channel(d_base, rx_mc, boxes)
    r_b = performance(d_base, h1, d1, hist_e10_base[beta20], tb)
    r_f = performance(d_filt, h1, d1, hist_e10_base[beta20], tb)
    r_c = performance(d_fec, h1, d1, hist_e10_base[beta20], tb)

    ha, da, _, _ = adr_channels(d_auto, rx_mc, boxes, include_refl=False)
    r_a = performance_adr(d_auto, ha, da, hist_e10_adr, tb)
    rw = wifi_evaluate(w_e10, rx_mc, HB, dyn, SF_MC, RNG, n_fad=1, use_spatial_emi=True)

    for k, r in (
        ("1. Base\nVLC", r_b),
        ("2. +Filtro\n+7 días", r_f),
        ("3. +RS\nFEC", r_c),
        ("4. +ADR+4×4\nAutónomo", r_a),
        ("5. Wi-Fi\n+FEC", rw),
    ):
        okD[k] += (r["ber"] <= BER_TH) & valid
        gpD[k].append(float(np.mean(r["goodput_mbps"][valid])))
        eeD[k].append(float(np.mean(r["ee_mbit_per_j"][valid])))
        if k in snr_store:
            snr_store[k].append(r["snr_db"][valid])

    hyb_ok = (r_a["ber"] <= BER_TH) | (rw["ber"] <= BER_TH)
    hyb_gp = np.maximum(r_a["goodput_mbps"][valid], rw["goodput_mbps"][valid])
    okD["6. Híbrido\nVLC+Wi-Fi"] += hyb_ok & valid
    gpD["6. Híbrido\nVLC+Wi-Fi"].append(float(np.mean(hyb_gp)))
    eeD["6. Híbrido\nVLC+Wi-Fi"].append(float(np.mean(hyb_gp) / (d_auto.P_comm_total_W + w_e10.P_comm_total_W)))

e10_av = {k: float(100.0 * np.mean((okD[k] / np.maximum(val_e10, 1.0)) >= 0.99)) for k in e10_keys}
e10_gp = {k: float(np.mean(gpD[k])) for k in e10_keys}
e10_ee = {k: float(np.mean(eeD[k])) for k in e10_keys}

fig, axs = plt.subplots(1, 3, figsize=(16.5, 4.9))
for k, c, lab in (
    ("1. Base\nVLC", "tab:red", "Base VLC (1-PD, 30 d)"),
    ("2. +Filtro\n+7 días", "tab:orange", "+Filtro azul + 7 d"),
    ("4. +ADR+4×4\nAutónomo", "tab:green", "VLC Autónomo (ADR + 4×4)"),
):
    s = np.sort(np.concatenate(snr_store[k]))
    axs[0].plot(s, np.arange(1, len(s) + 1) / len(s) * 100.0, c=c, lw=1.8, label=lab)
axs[0].axvline(SNR_TH_OOK_DB, color="k", ls="--", lw=1.0)
axs[0].axvline(SNR_TH_OOK_FEC_DB, color="tab:green", ls=":", lw=1.2)
axs[0].text(SNR_TH_OOK_DB + 0.6, 48, "13,5 dB\n(sin FEC)", fontsize=7.5)
axs[0].text(SNR_TH_OOK_FEC_DB - 9.0, 68, "10,0 dB\n(RS-FEC)", fontsize=7.5, color="darkgreen")
axs[0].set_xlim(-10, 45)
axs[0].set_xlabel("SNR intermedia [dB]")
axs[0].set_ylabel("CDF [%]")
axs[0].set_title("E10a · CDF de SNR en escenario combinado hostil")
axs[0].legend(fontsize=7.6, loc="upper left")

bar_colors = ["tab:red", "tab:orange", "gold", "tab:green", "tab:cyan", "tab:olive"]
x_k = np.arange(len(e10_keys))
vals_av = [e10_av[k] for k in e10_keys]
bars_av = axs[1].bar(x_k, vals_av, color=bar_colors, edgecolor="k", lw=0.6)
for b_, v in zip(bars_av, vals_av):
    axs[1].text(b_.get_x() + b_.get_width() / 2, v + 1.5, f"{v:.1f} %", ha="center", fontsize=7.8, fontweight="bold")
axs[1].set_xticks(x_k)
axs[1].set_xticklabels(e10_keys, fontsize=7.6)
axs[1].set_ylabel("% ubicaciones con disponibilidad ≥ 99 %")
axs[1].set_ylim(0, 112)
axs[1].set_title("E10b · Métrica 1 (Confiabilidad): Ablación VLC vs Wi-Fi")

vals_gp = [e10_gp[k] for k in e10_keys]
vals_ee = [e10_ee[k] for k in e10_keys]
bars_gp = axs[2].bar(x_k, vals_gp, color=bar_colors, edgecolor="k", lw=0.6)
for b_, v, ee_v in zip(bars_gp, vals_gp, vals_ee):
    axs[2].text(b_.get_x() + b_.get_width() / 2, v + 0.22, f"{v:.2f}\n({ee_v:.2f} M/J)", ha="center", fontsize=7.0)
axs[2].set_xticks(x_k)
axs[2].set_xticklabels(e10_keys, fontsize=7.6)
axs[2].set_ylabel("Goodput efectivo medio [Mb/s]")
axs[2].set_ylim(0, 11.6)
axs[2].set_title("E10c · Métrica 2 (Eficiencia): Goodput y Mbit/J")

fig.savefig(f"{OUT}/E10_diseno_mejorado.png")
plt.close(fig)

rows_e10 = []
for k in e10_keys:
    k_clean = k.replace("\n", " ")
    rows_e10.append([k_clean, f"{e10_av[k]:.1f}", f"{e10_gp[k]:.2f}", f"{e10_ee[k]:.3f}"])
    csv_rows.append(["E10", f"{k_clean} · Dispo ≥99 % [%] / Goodput [Mb/s] / EE [Mbit/J]",
                     f"{e10_av[k]:.1f} % / {e10_gp[k]:.2f} / {e10_ee[k]:.2f}"])
write_csv("E10_ablacion.csv", ["configuracion", "disponibilidad_99_%", "goodput_medio_Mbps", "eficiencia_Mbit_por_J"], rows_e10)
log("E10 listo")

# ===========================================================================
# Resumen final y registro de ejecución
# ===========================================================================
write_csv("resumen_resultados.csv", ["experimento", "indicador", "valor"], csv_rows)
for r_ in csv_rows:
    log(" | ".join(r_), timed=False)
log("Terminado")

with open(f"{OUT}/log_ejecucion.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(log_lines) + "\n")
