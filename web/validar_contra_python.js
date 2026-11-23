const V = require("./vlc_physics.js");
const ref = require("./referencia_python.json");
const P = V.defaults();
const grid = V.receiverGrid(P, 0.5);
const SB = V.STATIC_BOXES;
const free = []; for (let i = 0; i < grid.n; i++) free.push(!V.insideAny(grid.pts[3*i], grid.pts[3*i+1], SB));
const med = (a) => { const s = [...a].sort((x, y) => x - y); const n = s.length; return n % 2 ? s[(n-1)/2] : 0.5*(s[n/2-1]+s[n/2]); };
const pct = (a, q) => { const s = [...a].sort((x, y) => x - y); const pos = (s.length - 1) * q / 100, lo = Math.floor(pos); return s[lo] + (s[Math.min(lo+1, s.length-1)] - s[lo]) * (pos - lo); };
function run(name, p, boxes) {
  const t = Date.now();
  const brs = V.branches(p);
  const refl = V.reflections(p, grid, boxes, brs);
  const r = V.evaluateVLC(p, grid, boxes, refl);
  const idx = [...Array(grid.n).keys()].filter((i) => boxes.length ? free[i] : true);
  const snr = idx.map((i) => r.snr[i]), ber = idx.map((i) => r.ber[i]), gp = idx.map((i) => r.gp[i]);
  const js = { med: med(snr), p5: pct(snr, 5), cov: 100 * ber.filter((b) => b <= 1e-6).length / ber.length, gp: gp.reduce((a, b) => a + b, 0) / gp.length,
    snr_pts: [0, 200, 465, 700].map((i) => r.snr[i]) };
  const py = ref[name];
  console.log(`${name.padEnd(14)} med JS ${js.med.toFixed(2)} / PY ${py.med.toFixed(2)} | p5 ${js.p5.toFixed(2)} / ${py.p5.toFixed(2)} | cov ${js.cov.toFixed(1)} / ${py.cov.toFixed(1)} | gp ${js.gp.toFixed(2)} / ${py.gp.toFixed(2)} | pts ${js.snr_pts.map(v=>v.toFixed(1)).join(",")} vs ${py.snr_pts.map(v=>v.toFixed(1)).join(",")}  (${Date.now()-t} ms)`);
}
const R = (o) => Object.assign({}, P, o);
run("nominal", P, []);
run("hostil", R({ visibility: 10, days_dirty: 30, T_amb: 70, I_bg: 5100e-6 }), []);
run("estaticos", P, SB);
run("polvo5_fec", R({ visibility: 5, use_fec: true, Rb: 20e6, use_blue_filter: true }), SB);
run("adr4x4", Object.assign(V.optimizedGrid(P), { adr_enabled: true, visibility: 20, use_fec: true, T_amb: 70 }), SB);
run("led_crudo_30M", R({ f_3db_led: 3e6, Rb: 20e6 }), []);
// Wi-Fi
const w = Object.assign(V.wifiDefaults(), { use_fec: true });
const heavy = SB.filter((b) => b[4] >= 2);
const rg = V.rng(1); let ok = 0, n = 0, avs = 0, gps = 0;
for (let i = 0; i < grid.n; i++) { if (!free[i]) continue; const r = V.wifiPoint(w, grid.pts[3*i], grid.pts[3*i+1], 1, heavy, [], 0, rg, 400); n++; if (r.avail >= 0.99) ok++; avs += r.avail; gps += r.gp; }
console.log(`wifi           cov99 JS ${(100*ok/n).toFixed(1)} / PY ${ref.wifi.cov99.toFixed(1)} | av_mean ${(100*avs/n).toFixed(2)} / ${ref.wifi.av_mean.toFixed(2)} | gp ${(gps/n).toFixed(2)} / ${ref.wifi.gp.toFixed(2)}`);
console.log("emi pts JS", [0,200,465,700].map(i=>V.emiAt(w, grid.pts[3*i], grid.pts[3*i+1], 1).toFixed(2)).join(","), "PY", ref.wifi.emi_pts.map(v=>v.toFixed(2)).join(","));
const il = V.illuminance(P, grid); console.log(`illum mean JS ${il.mean.toFixed(1)} / PY ${ref.illum.mean.toFixed(1)} | U0 ${il.u0.toFixed(3)} / ${ref.illum.u0.toFixed(3)}`);
console.log("RS-FEC 1e-3 ->", V.rsPostFecBer(1e-3).toExponential(3), " 5e-4 ->", V.rsPostFecBer(5e-4).toExponential(3));
