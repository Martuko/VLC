import json, sys, pathlib, numpy as np
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
from dataclasses import replace
from vlc_core import *
P = Params(dA_wall=0.75)
SB = static_boxes()
rx, xs, ys = receiver_grid(P, 0.5)
free = ~inside_any(rx[:, :2], SB)
out = {}
def run(name, p, boxes, adr=False):
    beta = p.beta_ext
    if adr:
        ha, da, hl, tb = adr_channels(p, rx, boxes, (beta,), include_refl=True)
        r = performance_adr(p, ha, da, [h[beta] for h in hl], tb)
    else:
        h, d = los_channel(p, rx, boxes)
        hist, tb = reflection_histogram(p, rx, boxes, (beta,))
        r = performance(p, h, d, hist[beta], tb)
    area = free if boxes is not None else np.ones(len(rx), bool)
    out[name] = dict(med=float(np.median(r["snr_db"][area])), p5=float(np.percentile(r["snr_db"][area], 5)),
                     cov=float(100*np.mean(r["ber"][area] <= 1e-6)), gp=float(np.mean(r["goodput_mbps"][area])),
                     snr_pts=[float(r["snr_db"][i]) for i in (0, 200, 465, 700)])
run("nominal", P, None)
run("hostil", replace(P, visibility=10.0, days_dirty=30, T_amb=70, I_bg=5100e-6), None)
run("estaticos", P, SB)
run("polvo5_fec", replace(P, visibility=5.0, use_fec=True, Rb=20e6, use_blue_filter=True), SB)
run("adr4x4", replace(optimized_led_grid(P), adr_enabled=True, visibility=20.0, use_fec=True, T_amb=70), SB, adr=True)
run("led_crudo_30M", replace(P, f_3db_led=3e6, Rb=20e6), None)
w = WifiParams(use_fec=True)
rw = wifi_evaluate(w, rx, heavy_static_boxes(), None, np.zeros(len(rx)), np.random.default_rng(1), n_fad=400, use_spatial_emi=True)
av = np.mean(rw["ber"] <= 1e-6, axis=1)
out["wifi"] = dict(cov99=float(100*np.mean(av[free] >= 0.99)), av_mean=float(100*np.mean(av[free])), gp=float(np.mean(rw["goodput_mbps"][free])),
                   emi_pts=[float(x) for x in wifi_spatial_emi_dB(rx[[0, 200, 465, 700]])])
E = illuminance(P, rx, None, include_reflections=True)
out["illum"] = dict(mean=float(E.mean()), u0=float(E.min()/E.mean()))
out["pts"] = rx[[0, 200, 465, 700]].tolist()
json.dump(out, open(str(pathlib.Path(__file__).with_name("referencia_python.json")), "w"), indent=1)
print(json.dumps(out, indent=1))
