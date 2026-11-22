# -*- coding: utf-8 -*-
"""
demo_interactiva.py
===================
Demostración en vivo para la presentación oral. Abre una ventana con el mapa de
SNR de la nave y controles deslizantes para cambiar, en tiempo real:

  * visibilidad por polvo en suspensión [m]
  * días sin limpiar las ópticas
  * temperatura ambiente [°C]
  * corriente de luz ambiente I_bg [µA]
  * número de operarios en movimiento
  * tasa de bits [Mb/s]
  * campo de visión del receptor (FOV) [°]

Botón "Nuevos obstáculos": vuelve a sortear la posición de operarios y grúas.

Uso:  python demo_interactiva.py      (requiere numpy, scipy y matplotlib)
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.widgets import Slider, Button
from matplotlib.patches import Rectangle
from vlc_core import *

P0 = Params()
SB = static_boxes()
STEP = 0.5
rx, xs, ys = receiver_grid(P0, STEP)
free = ~inside_any(rx[:, :2], SB)
V_OPTS = np.array([5, 7, 10, 15, 20, 30, 50, 100, 200, 1e9])     # 1e9 ≈ aire limpio
FOV_OPTS = np.array([40, 50, 60, 70, 80])
rng = np.random.default_rng(7)

print("Precalculando reflexiones (unos segundos)...")
cache = {}
for fov in FOV_OPTS:
    p = replace(P0, fov_deg=float(fov))
    betas = tuple(3.912 / v for v in V_OPTS)
    hist, tb = reflection_histogram(p, rx, SB, betas)
    cache[fov] = hist
print("Listo.")

state = {"dyn": random_dynamic_boxes(rng, 10, 1, P0), "nw": 10}

fig = plt.figure(figsize=(11, 7.5))
ax = fig.add_axes([0.06, 0.36, 0.55, 0.6])
ax_txt = fig.add_axes([0.66, 0.36, 0.32, 0.6]); ax_txt.axis("off")
im = ax.imshow(np.zeros((len(ys), len(xs))), origin="lower", extent=[0, P0.L, 0, P0.W],
               cmap="RdYlGn", vmin=-10, vmax=50)
fig.colorbar(im, ax=ax, label="SNR [dB]  (13,5 dB → BER 10⁻⁶)")
ax.set_xlabel("x [m]"); ax.set_ylabel("y [m]")
patches = []

def mk(y, label, vmin, vmax, v0, step=None):
    a = fig.add_axes([0.18, y, 0.5, 0.025])
    return Slider(a, label, vmin, vmax, valinit=v0, valstep=step)

s_v = mk(0.29, "Visibilidad (índice)", 0, len(V_OPTS) - 1, len(V_OPTS) - 1, 1)
s_d = mk(0.255, "Días sin limpieza", 0, 120, 0, 1)
s_t = mk(0.22, "Temperatura [°C]", 0, 90, 25, 1)
s_b = mk(0.185, "Luz ambiente I_bg [µA]", 100, 6000, 1000, 100)
s_w = mk(0.15, "Operarios", 0, 40, 10, 1)
s_r = mk(0.115, "Tasa [Mb/s]", 1, 150, 50, 1)
s_f = mk(0.08, "FOV (índice: 40–80°)", 0, len(FOV_OPTS) - 1, 2, 1)
b_new = Button(fig.add_axes([0.75, 0.08, 0.18, 0.05]), "Nuevos obstáculos")

def update(_=None):
    V = V_OPTS[int(s_v.val)]
    fov = FOV_OPTS[int(s_f.val)]
    nw = int(s_w.val)
    if nw != state["nw"]:
        state["dyn"] = random_dynamic_boxes(rng, nw, nw // 10, P0); state["nw"] = nw
    p = replace(P0, visibility=(np.inf if V > 1e8 else float(V)), days_dirty=s_d.val, T_amb=s_t.val,
                I_bg=s_b.val * 1e-6, Rb=s_r.val * 1e6, fov_deg=float(fov))
    boxes = np.vstack([SB, state["dyn"]])
    h, d = los_channel(p, rx, boxes)
    hist = cache[fov][3.912 / V]
    r = performance(p, h, d, hist, tb)
    snr = np.where(free, r["snr_db"], np.nan)
    im.set_data(np.clip(snr, -10, 50).reshape(len(ys), len(xs)))
    for pt in patches:
        pt.remove()
    patches.clear()
    for b in SB:
        patches.append(ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2], fc="0.55", ec="k")))
    for b in state["dyn"]:
        c = "tab:orange" if (b[1] - b[0]) < 0.6 else "tab:red"
        patches.append(ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2], fc=c, ec="k", lw=0.5)))
    cov = 100 * np.mean(r["ber"][free] <= 1e-6)
    ax.set_title(f"Cobertura BER ≤ 10⁻⁶: {cov:.1f} %   ·   SNR mediana {np.nanmedian(snr):.1f} dB")
    txt = (f"Parámetros actuales\n\n"
           f"Visibilidad: {'aire limpio' if V > 1e8 else f'{V:.0f} m'}\n"
           f"β extinción: {p.beta_ext:.3f} 1/m\n"
           f"Transmitancia por suciedad: {p.T_soil:.2f}\n"
           f"Potencia LED relativa (T): {p.led_derating:.2f}\n"
           f"FOV: {fov}°  ·  Rb: {p.Rb/1e6:.0f} Mb/s\n"
           f"Operarios: {nw}  ·  grúas: {nw // 10}\n\n"
           f"Peor 5 % SNR: {np.nanpercentile(snr, 5):.1f} dB\n"
           f"BER mediana: {np.median(r['ber'][free]):.1e}")
    ax_txt.clear(); ax_txt.axis("off"); ax_txt.text(0, 1, txt, va="top", family="monospace", fontsize=10)
    fig.canvas.draw_idle()

def new_obs(_):
    state["dyn"] = random_dynamic_boxes(rng, state["nw"], state["nw"] // 10, P0); update()

for s in (s_v, s_d, s_t, s_b, s_w, s_r, s_f):
    s.on_changed(update)
b_new.on_clicked(new_obs)
update()
plt.show()
