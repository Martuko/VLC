/* validar_graficos.js — Prueba automatizada en Node de los 9 generadores dinámicos (E2–E11)
 * en web/vlc_charts.js sobre el motor físico web/vlc_physics.js. */
const assert = require("assert");
const V = require("./vlc_physics.js");
const C = require("./vlc_charts.js");

const p = Object.assign(V.defaults(), { visibility: 20, days_dirty: 30, T_amb: 70, I_bg: 5100e-6 });
const w = V.wifiDefaults();
const grid = V.receiverGrid(p, 0.5);
const dyn = V.randomDynamic(V.rng(7), 12, 3, p);
const boxes = V.STATIC_BOXES.concat(dyn);
const heavy = V.STATIC_BOXES.filter((b) => b[4] >= 2).concat(dyn.filter((b) => b[4] >= 2));
const workersB = dyn.filter((b) => b[4] < 2);
const freeIdx = [];
for (let i = 0; i < grid.n; i++) {
  if (!V.insideAny(grid.pts[3 * i], grid.pts[3 * i + 1], boxes)) freeIdx.push(i);
}
const sfRng = V.rng(2023), sfNorm = new Float64Array(grid.n).map(() => sfRng.normal());
const refl = V.reflections(p, grid, V.STATIC_BOXES, V.branches(p));
const losAll = V.losGrid(p, grid, boxes, V.branches(p));
const vlc = V.evaluateFromLos(p, grid, losAll, refl);
const rW = V.rng(18), wa = new Float32Array(grid.n), wg = new Float32Array(grid.n), we = new Float32Array(grid.n);
for (let i = 0; i < grid.n; i++) {
  const o = V.wifiPoint(w, grid.pts[3 * i], grid.pts[3 * i + 1], 1, heavy, workersB, sfNorm[i], rW, 40);
  wa[i] = o.avail; wg[i] = o.gp; we[i] = o.emi;
}
let cv = 0, gv = 0, cw = 0, gw = 0;
for (const i of freeIdx) {
  if (vlc.ber[i] <= 1e-6) cv++; gv += vlc.gp[i];
  if (wa[i] >= 0.99) cw++; gw += wg[i];
}
const last = { p, w, refl, losAll, boxes, heavy, workersB, vlc, wa, wg, we, freeIdx,
  covV: 100 * cv / freeIdx.length, covW: 100 * cw / freeIdx.length,
  gpV: gv / freeIdx.length, gpW: gw / freeIdx.length };

const reflMap = new Map(), losMap = new Map();
const evalConfigInstant = (pCfg, idxList) => {
  const kR = `${pCfg.led_xy.length}|${pCfg.adr_enabled}|${pCfg.adr_enabled ? pCfg.adr_fov_deg : pCfg.fov_deg}`;
  if (!reflMap.has(kR)) reflMap.set(kR, V.reflections(pCfg, grid, V.STATIC_BOXES, V.branches(pCfg)));
  if (!losMap.has(kR)) losMap.set(kR, V.losGrid(pCfg, grid, boxes, V.branches(pCfg)));
  return V.evaluateFromLos(pCfg, grid, losMap.get(kR), reflMap.get(kR), idxList);
};

const ctx = {
  V, p, w, grid, last, freeIdx, boxes, heavy, workersB, losAll, refl, sfNorm,
  seed: 7, workers: 12, agv: [6.25, 6.25], mcTab: null, evalConfigInstant,
};

const tabs = ["e2", "e3", "e4", "e5", "e6", "e7", "e9", "e11", "e8e10"];
for (const tab of tabs) {
  const t0 = Date.now();
  const data = C.computeTabData(tab, ctx);
  const dt = Date.now() - t0;
  assert(data.panels.length >= 1, `Tab ${tab} sin paneles`);
  for (const pn of data.panels) {
    if (pn.series) {
      for (const s of pn.series) {
        assert(s.xs.length > 0 && s.xs.length === s.ys.length, `${tab} longitud inválida`);
        for (const y of s.ys) assert(Number.isFinite(y), `${tab}/${s.label} contiene valor no finito: ${y}`);
      }
    }
    if (pn.bars) {
      for (const b of pn.bars) assert(Number.isFinite(b.value), `${tab} barra no finita`);
    }
  }
  console.log(`${tab.padEnd(6)} OK  (${data.panels.length} panel${data.panels.length > 1 ? "es" : ""}, ${dt} ms)  ${data.panels[0].title}`);
}

// Invariantes físicos específicos
const e4 = C.computeTabData("e4", ctx);
const covDust = e4.panels[1].series[0].ys;
for (let i = 1; i < covDust.length; i++) {
  assert(covDust[i] <= covDust[i - 1] + 1e-6, `E4 cobertura no monótona: ${covDust}`);
}

const e6 = C.computeTabData("e6", ctx);
const dSnrProt = e6.panels[0].series[0].ys;
const idx70 = e6.panels[0].series[0].xs.indexOf(70);
assert(dSnrProt[idx70] < -2.0, `E6 debe mostrar caída térmica real a 70°C, obtenido: ${dSnrProt[idx70]}`);

// E11 equilibrio con caso ideal (10 Mb/s VLC vs 9.31 Mb/s Wi-Fi, AP 12 W) -> ~1.56 W para 9 lum
const e11Clean = C.computeTabData("e11", Object.assign({}, ctx, {
  p: V.defaults(),
  w: Object.assign(V.wifiDefaults(), { P_ap_W: 12 }),
  last: Object.assign({}, last, { gpV: 10.0, gpW: 9.31 }),
}));
assert(e11Clean.panels[0].title.includes("1,56") || e11Clean.panels[0].title.includes("1,55"),
  `E11 equilibrio inesperado: ${e11Clean.panels[0].title}`);

console.log("Todas las pruebas de gráficos dinámicos E2–E11 pasaron correctamente.");
