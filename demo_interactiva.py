# -*- coding: utf-8 -*-
"""
demo_interactiva.py
===================
Demostración interactiva en vivo para la presentación oral del Proyecto 2.
Muestra en tiempo real el mapa de la nave industrial (15 × 15 × 6 m) y las dos
métricas oficiales (Confiabilidad vía BER y Eficiencia vía Goodput / Mbit/J),
permitiendo explorar:

  * Visibilidad por polvo en suspensión [m] (con extinción + scattering)
  * Días sin limpieza de ópticas (LDD + acoplamiento térmico e I_bg)
  * Temperatura ambiente [°C] (calculando temperatura de unión real Tj)
  * Luz ambiente I_bg [µA]
  * Operarios y grúas horquilla móviles (+ botón "Nuevos obstáculos")
  * Tasa de bits [Mb/s] (con ancho de banda realista del LED blanco)
  * Conmutadores en vivo:
      - Filtro óptico azul (25 MHz + rechazo solar/horno)
      - Codificación RS(255,239) FEC (IEEE 802.15.7)
      - Diseño mejorado ADR piramidal + grilla 4×4 EN 12464-1
      - Comparación en vivo: mapa de disponibilidad Wi-Fi 2,4 GHz bajo EMI de planta

Uso:
    python demo_interactiva.py
"""
from dataclasses import replace
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.widgets import Button, CheckButtons, Slider

from vlc_core import (
    EMI_SOURCES,
    Params,
    WifiParams,
    adr_channels,
    heavy_static_boxes,
    inside_any,
    los_channel,
    optimized_led_grid,
    performance,
    performance_adr,
    random_dynamic_boxes,
    receiver_grid,
    reflection_histogram,
    static_boxes,
    wifi_evaluate,
)


def build_app(show: bool = True):
    # En la demo interactiva se usa dA_wall = 0,75 m y paso 0,75 m para arranque rápido (~5–10 s)
    P0 = replace(Params(), dA_wall=0.75)
    P_OPT = replace(optimized_led_grid(P0), adr_enabled=True)
    SB = static_boxes()
    HB = heavy_static_boxes(min_height=2.0)
    STEP = 0.75
    rx, xs, ys = receiver_grid(P0, STEP)
    free_static = ~inside_any(rx[:, :2], SB)
    V_OPTS = [5.0, 7.0, 10.0, 15.0, 20.0, 50.0, np.inf]
    BETAS = tuple(sorted(set(0.0 if np.isinf(v) else 3.912 / v for v in V_OPTS)))
    rng = np.random.default_rng(7)
    sf_wifi = rng.standard_normal(len(rx))

    print("Precalculando reflexiones (base 1-PD y ADR 4×4, ~5–10 s)...", flush=True)
    hist_base, tb = reflection_histogram(P0, rx, SB, BETAS)
    _, _, hist_adr_list, _ = adr_channels(P_OPT, rx, SB, BETAS, include_refl=True)
    print("Listo.", flush=True)

    state = {
        "dyn": random_dynamic_boxes(rng, 10, 1, P0),
        "nw": 10,
        "blue": False,
        "fec": False,
        "adr_4x4": False,
        "wifi_mode": False,
    }

    fig = plt.figure(figsize=(12.2, 8.0))
    ax = fig.add_axes([0.05, 0.36, 0.54, 0.58])
    ax_txt = fig.add_axes([0.65, 0.44, 0.33, 0.50])
    ax_txt.axis("off")

    im = ax.imshow(
        np.zeros((len(ys), len(xs))),
        origin="lower",
        extent=[0, P0.L, 0, P0.W],
        cmap="RdYlGn",
        vmin=-5,
        vmax=45,
    )
    cbar = fig.colorbar(im, ax=ax, label="SNR intermedia [dB]  (13,5 dB → BER 10⁻⁶)")
    ax.set_xlabel("x [m]")
    ax.set_ylabel("y [m]")
    patches = []
    artists = []

    def mk_slider(y, label, vmin, vmax, v0, step=None):
        a = fig.add_axes([0.20, y, 0.45, 0.024])
        return Slider(a, label, vmin, vmax, valinit=v0, valstep=step)

    s_v = mk_slider(0.285, "Visibilidad (índice)", 0, len(V_OPTS) - 1, len(V_OPTS) - 1, 1)
    s_d = mk_slider(0.248, "Días sin limpieza", 0, 120, 0, 1)
    s_t = mk_slider(0.211, "Temperatura amb. [°C]", 0, 90, 25, 1)
    s_b = mk_slider(0.174, "Luz ambiente I_bg [µA]", 100, 6000, 1000, 100)
    s_w = mk_slider(0.137, "Operarios (+grúas)", 0, 40, 10, 1)
    s_r = mk_slider(0.100, "Tasa de bits [Mb/s]", 1, 50, 10, 1)

    ax_chk = fig.add_axes([0.68, 0.17, 0.28, 0.20])
    chk_labels = [
        "Filtro óptico azul (25 MHz)",
        "Codificación RS(255,239) FEC",
        "ADR piramidal + 4×4 EN12464",
        "Ver Wi-Fi 2,4 GHz (EMI planta)",
    ]
    chk = CheckButtons(ax_chk, chk_labels, [False, False, False, False])
    b_new = Button(fig.add_axes([0.68, 0.095, 0.28, 0.05]), "Sortear nuevos obstáculos")

    def update(_=None):
        V = V_OPTS[int(s_v.val)]
        beta = 0.0 if np.isinf(V) else 3.912 / V
        nw = int(s_w.val)
        if nw != state["nw"]:
            state["dyn"] = random_dynamic_boxes(rng, nw, nw // 10, P0)
            state["nw"] = nw

        boxes = np.vstack([SB, state["dyn"]]) if len(state["dyn"]) else SB
        free_now = ~inside_any(rx[:, :2], boxes)

        for pt in patches:
            pt.remove()
        patches.clear()
        for ar in artists:
            ar.remove()
        artists.clear()

        for b in SB:
            patches.append(ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2], fc="0.55", ec="k", lw=0.7)))
        for b in state["dyn"]:
            c = "tab:orange" if (b[1] - b[0]) < 0.7 else "tab:red"
            patches.append(ax.add_patch(Rectangle((b[0], b[2]), b[1] - b[0], b[3] - b[2], fc=c, ec="k", lw=0.5)))

        if state["wifi_mode"]:
            w_cfg = WifiParams(Rb_target=float(s_r.val) * 1e6, use_fec=state["fec"])
            rw = wifi_evaluate(w_cfg, rx, HB, state["dyn"], sf_wifi, np.random.default_rng(11),
                               n_fad=80, use_spatial_emi=True)
            av_pct = 100.0 * np.mean(rw["ber"] <= 1e-6, axis=1)
            data_map = np.where(free_now, av_pct, np.nan)
            im.set_clim(40, 100)
            im.set_data(np.clip(data_map, 40, 100).reshape(len(ys), len(xs)))
            cbar.set_label("Disponibilidad Wi-Fi (BER ≤ 10⁻⁶) [%]")

            xs_e = [v[0] for v in EMI_SOURCES.values()]
            ys_e = [v[1] for v in EMI_SOURCES.values()]
            artists.extend(ax.plot(xs_e, ys_e, "o", ms=7.5, mfc="none", mec="crimson", mew=1.8))
            artists.extend(ax.plot(w_cfg.ap_pos[0], w_cfg.ap_pos[1], "^", ms=9, mfc="tab:cyan", mec="k"))

            cov99 = 100.0 * np.mean(av_pct[free_now] >= 99.0)
            gp_m = float(np.mean(rw["goodput_mbps"][free_now]))
            ee_m = gp_m / w_cfg.P_comm_total_W
            ax.set_title(f"Wi-Fi 2,4 GHz (EMI espacial) · Dispo ≥99 %: {cov99:.1f} % · Gp: {gp_m:.2f} Mb/s")
            txt = (
                "MODO: Wi-Fi 2,4 GHz (EMI Planta)\n"
                "--------------------------------\n"
                f"Tasa objetivo : {w_cfg.Rb_target/1e6:.0f} Mb/s\n"
                f"RS(255,239)   : {'SÍ' if w_cfg.use_fec else 'NO'}\n"
                f"EMI espacial  : {rw['emi_dB'].min():.0f}–{rw['emi_dB'].max():.0f} dB\n"
                f"Operarios     : {nw} (grúas: {nw // 10})\n\n"
                "MÉTRICAS OFICIALES\n"
                f"1. Confiab. (≥99%): {cov99:.1f} %\n"
                f"2. Goodput medio  : {gp_m:.2f} Mb/s\n"
                f"   Eficiencia     : {ee_m:.2f} Mbit/J\n\n"
                "Círculos rojos = fuentes EMI\n"
                "(horno, soldadura, VFD)"
            )
        else:
            base_cfg = P_OPT if state["adr_4x4"] else P0
            p = replace(
                base_cfg,
                visibility=V,
                days_dirty=float(s_d.val),
                T_amb=float(s_t.val),
                I_bg=float(s_b.val) * 1e-6,
                Rb=float(s_r.val) * 1e6,
                use_blue_filter=state["blue"],
                use_fec=state["fec"],
            )
            if state["adr_4x4"]:
                ha, da, _, _ = adr_channels(p, rx, boxes, include_refl=False)
                hi_b = [hb[beta] for hb in hist_adr_list]
                r = performance_adr(p, ha, da, hi_b, tb)
            else:
                h, d = los_channel(p, rx, boxes)
                r = performance(p, h, d, hist_base[beta], tb)

            snr = np.where(free_now, r["snr_db"], np.nan)
            im.set_clim(-5, 45)
            im.set_data(np.clip(snr, -5, 45).reshape(len(ys), len(xs)))
            cbar.set_label("SNR intermedia [dB]  (13,5 dB sin FEC / 10,0 dB RS-FEC)")
            artists.extend(ax.plot(p.leds[:, 0], p.leds[:, 1], "*", ms=8, mfc="gold", mec="k"))

            cov = 100.0 * np.mean(r["ber"][free_now] <= 1e-6)
            gp_m = float(np.mean(r["goodput_mbps"][free_now]))
            ee_m = gp_m / p.P_comm_total_W
            ax.set_title(f"VLC {'ADR+4×4' if state['adr_4x4'] else 'Base 3×3'} · Cobertura: {cov:.1f} % · Goodput: {gp_m:.2f} Mb/s")
            txt = (
                f"MODO: VLC ({'ADR + 4×4' if state['adr_4x4'] else '1-PD + 3×3'})\n"
                "--------------------------------\n"
                f"Visibilidad   : {'aire limpio' if np.isinf(V) else f'{V:.0f} m'}\n"
                f"β extinción   : {p.beta_ext:.3f} 1/m (c/scatt.)\n"
                f"Suciedad T_s  : {p.T_soil:.2f} ({s_d.val:.0f} d)\n"
                f"Temp. unión Tj: {p.T_j:.1f} °C (der={p.led_derating:.2f})\n"
                f"BW LED efec.  : {p.f_3db_effective/1e6:.0f} MHz | Rb={p.Rb/1e6:.0f} Mb/s\n"
                f"RS(255,239)   : {'SÍ' if p.use_fec else 'NO'}\n"
                f"Operarios     : {nw} (grúas: {nw // 10})\n\n"
                "MÉTRICAS OFICIALES\n"
                f"1. Cobertura BER≤10⁻⁶: {cov:.1f} %\n"
                f"   BER mediana       : {np.median(r['ber'][free_now]):.1e}\n"
                f"2. Goodput medio     : {gp_m:.2f} Mb/s\n"
                f"   Eficiencia        : {ee_m:.2f} Mbit/J\n\n"
                f"SNR mediana (diag.)  : {np.nanmedian(snr):.1f} dB"
            )

        ax_txt.clear()
        ax_txt.axis("off")
        ax_txt.text(0.0, 1.0, txt, va="top", family="monospace", fontsize=9.0)
        fig.canvas.draw_idle()

    def on_check(label):
        if label == chk_labels[0]:
            state["blue"] = not state["blue"]
        elif label == chk_labels[1]:
            state["fec"] = not state["fec"]
        elif label == chk_labels[2]:
            state["adr_4x4"] = not state["adr_4x4"]
        elif label == chk_labels[3]:
            state["wifi_mode"] = not state["wifi_mode"]
        update()

    def on_new_obs(_):
        state["dyn"] = random_dynamic_boxes(rng, state["nw"], state["nw"] // 10, P0)
        update()

    for s in (s_v, s_d, s_t, s_b, s_w, s_r):
        s.on_changed(update)
    chk.on_clicked(on_check)
    b_new.on_clicked(on_new_obs)

    update()
    if show:
        plt.show()
    return fig, update, on_check, chk_labels


if __name__ == "__main__":
    build_app(show=True)
