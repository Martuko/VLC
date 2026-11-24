/* vlc_charts.js — Motor de gráficos dinámicos en Canvas 2D (sin dependencias externas)
 * para las pestañas E2–E11 del simulador 3D VLC industrial.
 * Se inyecta en docs/index.html mediante web/build.py y es importable en Node para validación. */
const VLCCharts = (() => {
  "use strict";

  const VIS_SWEEP = [200, 100, 50, 30, 20, 15, 10, 7, 5];
  const DAYS_SWEEP = [0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100, 110, 120];
  const TEMP_SWEEP = [0, 10, 20, 25, 35, 45, 55, 65, 70, 80, 90];
  const RATE_SWEEP = [1, 2, 5, 10, 15, 20, 30, 50, 75, 100];
  const EMI_SWEEP = [0, 10, 15, 20, 25, 30, 35, 40];

  const fmt = (v, d = 1) => (Number.isFinite(v) ? v.toFixed(d).replace(".", ",") : "–");
  function sciText(v) {
    if (!(v > 0)) return "0";
    const e = Math.floor(Math.log10(v)), m = v / Math.pow(10, e);
    return `${fmt(m, 1)}e${e}`;
  }
  function median(arr) {
    if (!arr.length) return NaN;
    const s = Array.from(arr).sort((a, b) => a - b);
    const n = s.length;
    return n % 2 ? s[(n - 1) >> 1] : 0.5 * (s[n / 2 - 1] + s[n / 2]);
  }
  function percentile(arr, q) {
    if (!arr.length) return NaN;
    const s = Array.from(arr).sort((a, b) => a - b);
    const pos = (s.length - 1) * q / 100, lo = Math.floor(pos), hi = Math.min(lo + 1, s.length - 1);
    return s[lo] + (s[hi] - s[lo]) * (pos - lo);
  }
  function sampleIndices(freeIdx, maxN = 180) {
    if (freeIdx.length <= maxN) return freeIdx;
    const step = freeIdx.length / maxN, out = [];
    for (let k = 0; k < maxN; k++) out.push(freeIdx[Math.floor(k * step)]);
    return out;
  }
  function statsVLC(res) {
    const n = res.ber.length;
    if (!n) return { cov: 0, gp: 0, snrMed: NaN, logBerMed: 0, logBerP95: 0 };
    let ok = 0, gpSum = 0;
    const lb = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      if (res.ber[i] <= 1e-6) ok++;
      gpSum += res.gp[i];
      lb[i] = Math.log10(Math.max(res.ber[i], 1e-30));
    }
    return {
      cov: 100 * ok / n,
      gp: gpSum / n,
      snrMed: median(res.snr),
      logBerMed: median(lb),
      logBerP95: percentile(lb, 95),
    };
  }

  // =========================================================================
  // 1. GENERADORES PUROS DE DATOS POR PESTAÑA (testeables sin DOM en Node)
  // =========================================================================
  function computeTabData(tabId, ctx) {
    const V = ctx.V, p = ctx.p, w = ctx.w, grid = ctx.grid;
    const freeIdx = ctx.freeIdx, idxS = sampleIndices(freeIdx, 180);
    const losAll = ctx.losAll, refl = ctx.refl;
    const P = (o) => Object.assign({}, p, o);

    // --- E2: Respuesta impulsiva en la posición actual del AGV ---
    if (tabId === "e2") {
      const [ax, ay] = ctx.agv;
      const brs = V.branches(p);
      const one = { pts: new Float64Array([ax, ay, p.z_rx]), n: 1 };
      const rOne = ctx.probeRefl || V.reflections(p, one, V.STATIC_BOXES, brs);
      const losOne = ctx.probeLos || V.losChannel(p, ax, ay, p.z_rx, ctx.boxes, brs);
      const imp = V.pointImpulse(p, brs, losOne, rOne, 0);
      const tNs = [], hVals = [];
      let maxH = 0;
      for (let k = 0; k < imp.hr.length; k++) {
        const t = k * imp.dtNs;
        if (t >= 10 && t <= 90) {
          tNs.push(t);
          hVals.push(imp.hr[k]);
          if (imp.hr[k] > maxH) maxH = imp.hr[k];
        }
      }
      const tauTxt = Number.isFinite(imp.tauRmsNs) ? `τ_RMS = ${fmt(imp.tauRmsNs, 1)} ns` : "sin señal";
      return {
        tabId,
        panels: [{
          kind: "impulse",
          title: `E2 · Respuesta impulsiva en el AGV (${fmt(ax, 1)}, ${fmt(ay, 1)}) m · ${imp.losArrivals.length} lum. visibles · ${tauTxt}`,
          xLabel: "Tiempo [ns]",
          yLabel: "h(t) reflexiones [1/s]",
          xMin: 10, xMax: 90,
          tNs, hVals,
          losArrivals: imp.losArrivals,
          t0Ns: imp.t0Ns,
          maxH: Math.max(maxH, 1e-6),
        }],
      };
    }

    // --- E3: CDF de BER en el área libre actual ---
    if (tabId === "e3") {
      const lbCur = [], lbIdeal = [];
      for (const i of freeIdx) lbCur.push(Math.log10(Math.max(ctx.last.vlc.ber[i], 1e-30)));
      lbCur.sort((a, b) => a - b);
      // Referencia ideal rápida (sin polvo/suciedad/calor/sol, con la geometría actual de luminarias)
      const pIdeal = P({ visibility: Infinity, days_dirty: 0, T_amb: 25, I_bg: 1000e-6 });
      const rIdeal = V.evaluateFromLos(pIdeal, grid, losAll, refl, idxS);
      for (let i = 0; i < rIdeal.ber.length; i++) lbIdeal.push(Math.log10(Math.max(rIdeal.ber[i], 1e-30)));
      lbIdeal.sort((a, b) => a - b);
      const toCdf = (sorted) => {
        const xs = [], ys = [], step = Math.max(1, Math.floor(sorted.length / 120));
        for (let i = 0; i < sorted.length; i += step) { xs.push(sorted[i]); ys.push(100 * (i + 1) / sorted.length); }
        if (sorted.length && xs[xs.length - 1] !== sorted[sorted.length - 1]) {
          xs.push(sorted[sorted.length - 1]); ys.push(100);
        }
        return { xs, ys };
      };
      const cCur = toCdf(lbCur), cId = toCdf(lbIdeal);
      return {
        tabId,
        panels: [{
          kind: "xy",
          title: `E3 · Distribución acumulada (CDF) de BER · Cobertura actual: ${fmt(ctx.last.covV, 1)} %`,
          xLabel: "log10(BER)",
          yLabel: "% área libre con BER ≤ x",
          xScale: "linear", xMin: -30, xMax: 0, yMin: 0, yMax: 103,
          vLines: [{ x: -6, label: "BER = 10⁻⁶", color: "fail", dashed: true }],
          series: [
            { label: "Configuración actual", xs: cCur.xs, ys: cCur.ys, color: "ok", width: 2.4 },
            { label: "Referencia limpia (25 °C, 0 d)", xs: cId.xs, ys: cId.ys, color: "muted", width: 1.6, dashed: true },
          ],
        }],
      };
    }

    // --- E4: Polvo en suspensión (BER + Cobertura vs Visibilidad) ---
    if (tabId === "e4") {
      const curRate = Math.round(p.Rb / 1e6);
      const rates = [curRate];
      for (const cand of [10, 30, 2]) if (rates.length < 3 && !rates.includes(cand)) rates.push(cand);
      const palette = ["ok", "blue", "hivis"];
      const sBer = [], sCov = [];
      rates.forEach((rbM, k) => {
        const med = [], p95 = [], cov = [];
        for (const vis of VIS_SWEEP) {
          const st = statsVLC(V.evaluateFromLos(P({ visibility: vis, Rb: rbM * 1e6 }), grid, losAll, refl, idxS));
          med.push(st.logBerMed); p95.push(st.logBerP95); cov.push(st.cov);
        }
        const tag = k === 0 ? `${rbM} Mb/s (actual)` : `${rbM} Mb/s`;
        sBer.push({ label: `${tag} · med.`, xs: VIS_SWEEP, ys: med, color: palette[k], width: k === 0 ? 2.4 : 1.7, dots: true });
        if (k === 0) sBer.push({ label: `${tag} · peor 5 %`, xs: VIS_SWEEP, ys: p95, color: palette[k], width: 1.5, dashed: true, dots: true });
        sCov.push({ label: tag, xs: VIS_SWEEP, ys: cov, color: palette[k], width: k === 0 ? 2.4 : 1.7, dots: true });
      });
      const marker = Number.isFinite(p.visibility)
        ? { x: p.visibility, label: `Actual: V=${p.visibility} m` }
        : { note: "Actual: aire limpio (V = ∞)" };
      return {
        tabId,
        panels: [
          {
            kind: "xy",
            title: "E4a · BER vs concentración de polvo (extinción + scattering)",
            xLabel: "Visibilidad V [m] (menor = más polvo)", yLabel: "log10(BER)",
            xScale: "log-inv", xTicks: [200, 100, 50, 20, 10, 5], yMin: -30, yMax: 0,
            hLines: [{ y: -6, label: "10⁻⁶", color: "fail", dashed: true }],
            currentMarker: marker, series: sBer,
          },
          {
            kind: "xy",
            title: "E4b · Cobertura confiable vs concentración de polvo",
            xLabel: "Visibilidad V [m] (menor = más polvo)", yLabel: "Cobertura BER ≤ 10⁻⁶ [%]",
            xScale: "log-inv", xTicks: [200, 100, 50, 20, 10, 5], yMin: 0, yMax: 104,
            currentMarker: marker, series: sCov,
          },
        ],
      };
    }

    // --- E5: Suciedad acumulada en ópticas ---
    if (tabId === "e5") {
      const visCases = [{ v: p.visibility, label: Number.isFinite(p.visibility) ? `Actual (V=${p.visibility} m)` : "Actual (aire limpio)", color: "ok", width: 2.5 }];
      if (Number.isFinite(p.visibility)) visCases.push({ v: Infinity, label: "Aire limpio", color: "blue", width: 1.6, dashed: true });
      if (p.visibility !== 20) visCases.push({ v: 20, label: "Polvo V = 20 m", color: "safety", width: 1.6 });
      if (p.visibility !== 10) visCases.push({ v: 10, label: "Polvo V = 10 m", color: "fail", width: 1.6 });
      const series = visCases.slice(0, 3).map((vc) => {
        const ys = DAYS_SWEEP.map((d) => statsVLC(V.evaluateFromLos(P({ visibility: vc.v, days_dirty: d }), grid, losAll, refl, idxS)).cov);
        return { label: vc.label, xs: DAYS_SWEEP, ys, color: vc.color, width: vc.width, dashed: !!vc.dashed, dots: true };
      });
      return {
        tabId,
        panels: [{
          kind: "xy",
          title: `E5 · Efecto de suciedad acumulada LDD + acople térmico/I_bg (${Math.round(p.Rb / 1e6)} Mb/s)`,
          xLabel: "Días sin limpiar ópticas", yLabel: "Cobertura BER ≤ 10⁻⁶ [%]",
          xScale: "linear", xMin: 0, xMax: 120, yMin: 0, yMax: 104,
          currentMarker: { x: p.days_dirty, label: `Actual: ${p.days_dirty} d` },
          series,
        }],
      };
    }

    // --- E6: Temperatura (ΔSNR + descomposición de ruido en el AGV) ---
    if (tabId === "e6") {
      const refProt = statsVLC(V.evaluateFromLos(P({ T_amb: 25 }), grid, losAll, refl, idxS)).snrMed;
      const refNoD = statsVLC(V.evaluateFromLos(P({ T_amb: 25, T_j_derate: null }), grid, losAll, refl, idxS)).snrMed;
      const yProt = [], yNoD = [];
      for (const T of TEMP_SWEEP) {
        yProt.push(statsVLC(V.evaluateFromLos(P({ T_amb: T }), grid, losAll, refl, idxS)).snrMed - refProt);
        yNoD.push(statsVLC(V.evaluateFromLos(P({ T_amb: T, T_j_derate: null }), grid, losAll, refl, idxS)).snrMed - refNoD);
      }
      const minY = Math.min(-1, ...yProt, ...yNoD);
      const tDerateAmb = 85 - p.delta_T_j_nom * (1 + p.k_soil_th * (1 - V.D.Tstx(p)));

      // E6b: descomposición de ruido en la posición actual del AGV
      const [ax, ay] = ctx.agv;
      const brs = V.branches(p);
      const one = { pts: new Float64Array([ax, ay, p.z_rx]), n: 1 };
      const rOne = ctx.probeRefl || V.reflections(p, one, V.STATIC_BOXES, brs);
      const losOne = ctx.probeLos || V.losChannel(p, ax, ay, p.z_rx, ctx.boxes, brs);
      const tempsBar = [25, 50, 70, 85];
      const shotVals = [], thVals = [], dkVals = [], barLabels = [];
      let activeGroup = 0, bestDiff = Infinity;
      tempsBar.forEach((T, idx) => {
        if (Math.abs(T - p.T_amb) < bestDiff) { bestDiff = Math.abs(T - p.T_amb); activeGroup = idx; }
        const pt = P({ T_amb: T });
        const e = V.eyeTerms(pt, brs[0], losOne[0].h, losOne[0].d, rOne ? rOne.hist[0] : null, 0, rOne ? rOne.dt : 0, rOne ? rOne.nb : 0);
        const nb = V.noiseBreakdown(pt, e.Pdc, brs[0].A);
        shotVals.push(Math.max(nb.sShotSignalBg, 1e-22));
        thVals.push(Math.max(nb.sTh, 1e-22));
        dkVals.push(Math.max(nb.sDark, 1e-22));
        barLabels.push(`${T}°C\nTj=${Math.round(V.D.Tj(pt))}°`);
      });
      return {
        tabId,
        panels: [
          {
            kind: "xy",
            title: "E6a · Pérdida térmica de SNR mediana (respecto a 25 °C)",
            xLabel: "Temperatura ambiente T_amb [°C]", yLabel: "ΔSNR mediana [dB]",
            xScale: "linear", xMin: 0, xMax: 90, yMin: Math.floor(minY * 1.15 - 0.3), yMax: Math.max(0.8, ...yProt) + 0.3,
            vLines: [{ x: tDerateAmb, label: "Tj=85°C", color: "muted", dashed: true }],
            currentMarker: { x: p.T_amb, label: `Actual: ${p.T_amb} °C` },
            series: [
              { label: "Con derating driver (Tj > 85 °C)", xs: TEMP_SWEEP, ys: yProt, color: "ok", width: 2.4, dots: true },
              { label: "Sin derating (solo droop LED)", xs: TEMP_SWEEP, ys: yNoD, color: "blue", width: 1.7, dashed: true, dots: true },
            ],
          },
          {
            kind: "bars-grouped-log",
            title: `E6b · Componentes de ruido en el AGV (${fmt(ax, 1)}, ${fmt(ay, 1)}) m`,
            yLabel: "Varianza de ruido [A²]",
            yMinExp: -20, yMaxExp: -13,
            labels: barLabels,
            activeGroup,
            groups: [
              { label: "Shot (señal+amb.)", values: shotVals, color: "blue" },
              { label: "Térmico", values: thVals, color: "hivis" },
              { label: "Corriente oscura", values: dkVals, color: "ok" },
            ],
          },
        ],
      };
    }

    // --- E7: Tasa de bits, BW del LED y Goodput ---
    if (tabId === "e7") {
      const cfgs = [
        { label: "Configuración actual", ov: {}, color: "ok", width: 2.5, dashed: false },
        { label: "LED sin ecualizar (3 MHz)", ov: { f_3db_led: 3e6, use_blue_filter: false }, color: "muted", width: 1.6, dashed: true },
        { label: "Filtro azul (25 MHz) + RS-FEC", ov: { use_blue_filter: true, use_fec: true }, color: "purple", width: 1.9, dashed: false },
      ];
      const sCov = [], sGp = [];
      for (const c of cfgs) {
        const covs = [], gps = [];
        for (const rbM of RATE_SWEEP) {
          const st = statsVLC(V.evaluateFromLos(P(Object.assign({}, c.ov, { Rb: rbM * 1e6 })), grid, losAll, refl, idxS));
          covs.push(st.cov); gps.push(st.gp);
        }
        sCov.push({ label: c.label, xs: RATE_SWEEP, ys: covs, color: c.color, width: c.width, dashed: c.dashed, dots: true });
        sGp.push({ label: c.label, xs: RATE_SWEEP, ys: gps, color: c.color, width: c.width, dashed: c.dashed, dots: true });
      }
      const maxGp = Math.max(12, ...sGp.flatMap((s) => s.ys));
      const curRb = p.Rb / 1e6;
      return {
        tabId,
        panels: [
          {
            kind: "xy",
            title: "E7a · Confiabilidad vs tasa de bits (ancho de banda LED e ISI)",
            xLabel: "Tasa de bits OOK [Mb/s]", yLabel: "Cobertura BER ≤ 10⁻⁶ [%]",
            xScale: "log", xTicks: [1, 2, 5, 10, 20, 50, 100], yMin: 0, yMax: 104,
            currentMarker: { x: curRb, label: `Actual: ${fmt(curRb, 0)} Mb/s` },
            series: sCov,
          },
          {
            kind: "xy",
            title: "E7b · Eficiencia (Goodput útil) vs tasa de transmisión",
            xLabel: "Tasa de bits OOK [Mb/s]", yLabel: "Goodput medio [Mb/s]",
            xScale: "log", xTicks: [1, 2, 5, 10, 20, 50, 100], yMin: 0, yMax: Math.ceil(maxGp * 1.12),
            currentMarker: { x: curRb, label: `Actual: ${fmt(curRb, 0)} Mb/s` },
            series: sGp,
          },
        ],
      };
    }

    // --- E9: Comparación VLC vs Wi-Fi (EMI impulsiva y Polvo) ---
    if (tabId === "e9") {
      const idxW = sampleIndices(freeIdx, 120);
      const sw = V.wifiSweepEmi(w, grid, idxW, ctx.heavy, ctx.workersB, ctx.sfNorm, EMI_SWEEP, 31 + ctx.seed, 30);
      const vlcCurCov = ctx.last.covV;
      const vlcEmiLine = EMI_SWEEP.map(() => vlcCurCov);
      const emiMed = median(idxW.map((i) => ctx.last.we[i]));

      const vlcDust = VIS_SWEEP.map((vis) => statsVLC(V.evaluateFromLos(P({ visibility: vis }), grid, losAll, refl, idxS)).cov);
      const wifiCurLine = VIS_SWEEP.map(() => ctx.last.covW);
      const dustMarker = Number.isFinite(p.visibility)
        ? { x: p.visibility, label: `Actual: V=${p.visibility} m` }
        : { note: "Actual: aire limpio" };

      return {
        tabId,
        panels: [
          {
            kind: "xy",
            title: `E9a · Interferencia EMI Middleton Class A (${Math.round(w.Rb_target / 1e6)} Mb/s)`,
            xLabel: "EMI impulsiva sobre ruido térmico [dB]", yLabel: "Confiabilidad [%]",
            xScale: "linear", xMin: 0, xMax: 40, yMin: 0, yMax: 105,
            currentMarker: { x: Math.min(40, Math.max(0, emiMed)), label: `EMI nave ~${fmt(emiMed, 0)} dB` },
            annotations: [
              { x: 4, y: 5, text: "Oficina" },
              { x: 14, y: 5, text: "Planta" },
              { x: 24, y: 5, text: "VFD" },
              { x: 35, y: 5, text: "Horno/sold." },
            ],
            series: [
              { label: "VLC actual (inmune EMI)", xs: EMI_SWEEP, ys: vlcEmiLine, color: "ok", width: 2.5, dots: true },
              { label: "Wi-Fi (+RS-FEC)", xs: EMI_SWEEP, ys: sw.covFec, color: "rf", width: 2.1, dots: true },
              { label: "Wi-Fi (sin FEC)", xs: EMI_SWEEP, ys: sw.covRaw, color: "rf", width: 1.5, dashed: true, dots: true },
            ],
          },
          {
            kind: "xy",
            title: "E9b · Polvo en suspensión (complementariedad VLC vs Wi-Fi)",
            xLabel: "Visibilidad V [m] (menor = más polvo)", yLabel: "Confiabilidad [%]",
            xScale: "log-inv", xTicks: [200, 100, 50, 20, 10, 5], yMin: 0, yMax: 105,
            currentMarker: dustMarker,
            series: [
              { label: "VLC actual", xs: VIS_SWEEP, ys: vlcDust, color: "ok", width: 2.5, dots: true },
              { label: `Wi-Fi (EMI nave, ${w.use_fec ? "+FEC" : "sin FEC"})`, xs: VIS_SWEEP, ys: wifiCurLine, color: "rf", width: 2.1, dashed: true, dots: true },
            ],
          },
        ],
      };
    }

    // --- E11: Sensibilidad de la eficiencia energética (Mbit/J) ---
    if (tabId === "e11") {
      const xsW = [];
      for (let x = 0.15; x <= 3.001; x += 0.05) xsW.push(+x.toFixed(2));
      const gpV = ctx.last.gpV, gpW = ctx.last.gpW;
      const prx1 = p.P_rx_vlc, prxAdr = p.P_rx_vlc * 1.6;
      const y9 = xsW.map((pw) => gpV / (9 * pw + prx1));
      const y16 = xsW.map((pw) => gpV / (16 * pw + prxAdr));
      const eeApCur = gpW / (w.P_ap_W + w.P_rx_wifi_W);
      const eeAp6 = gpW / (6 + w.P_rx_wifi_W);
      const eeAp15 = gpW / (15 + w.P_rx_wifi_W);

      const nl = p.led_xy.length, prxCur = p.P_rx_vlc * (p.adr_enabled ? 1.6 : 1);
      const beCur = gpW > 0 ? (gpV * (w.P_ap_W + w.P_rx_wifi_W) / gpW - prxCur) / nl : NaN;
      const beTxt = Number.isFinite(beCur) && beCur > 0
        ? `Equilibrio (${nl} lum. vs AP ${fmt(w.P_ap_W, 0)} W): ${fmt(beCur, 2)} W/lum.`
        : "Sin equilibrio positivo";

      return {
        tabId,
        panels: [{
          kind: "xy",
          title: `E11 · Sensibilidad de eficiencia energética · ${beTxt}`,
          xLabel: "Potencia eléctrica extra por luminaria para modular [W]",
          yLabel: "Eficiencia incremental [Mbit/J]",
          xScale: "linear", xMin: 0, xMax: 3.0, yMin: 0, yMax: 3.0,
          currentMarker: { x: p.P_comm_extra_per_led, label: `Actual: ${fmt(p.P_comm_extra_per_led, 2)} W` },
          hLines: [
            { y: eeApCur, label: `Wi-Fi AP ${fmt(w.P_ap_W, 0)} W (${fmt(eeApCur, 2)})`, color: "rf", dashed: false },
            { y: eeAp6, label: `AP 6 W (${fmt(eeAp6, 2)})`, color: "muted", dashed: true },
            { y: eeAp15, label: `AP 15 W (${fmt(eeAp15, 2)})`, color: "muted", dashed: true },
          ],
          series: [
            { label: `VLC 9 lum. 1-PD (${fmt(gpV, 2)} Mb/s)`, xs: xsW, ys: y9, color: "blue", width: p.led_xy.length === 9 ? 2.6 : 1.7 },
            { label: `VLC 16 lum. ADR (${fmt(gpV, 2)} Mb/s)`, xs: xsW, ys: y16, color: "ok", width: p.led_xy.length === 16 ? 2.6 : 1.7 },
          ],
        }],
      };
    }

    // --- E8/E10: Ablación instantánea o Monte Carlo bajo demanda ---
    if (tabId === "e8e10") {
      if (ctx.mcTab && ctx.mcTab.valid) {
        const mc = ctx.mcTab;
        return {
          tabId,
          mcActive: true,
          panels: [
            {
              kind: "xy",
              title: `E8 · Disponibilidad ≥ 99 % vs tránsito (${mc.nRuns} sorteos Monte Carlo)`,
              xLabel: "Operarios en la nave (+1 grúa cada 10)", yLabel: "Ubicaciones con dispo. ≥ 99 % [%]",
              xScale: "linear", xMin: 0, xMax: 40, yMin: 0, yMax: 105,
              currentMarker: { x: ctx.workers, label: `Actual: ${ctx.workers} op.` },
              series: mc.e8Series,
            },
            {
              kind: "bars",
              title: `E10 · Ablación Monte Carlo (dispo. ≥ 99 %, ${mc.nw10 ?? ctx.workers} op.)`,
              yLabel: "Ubicaciones ≥ 99 % [%]",
              yMin: 0, yMax: 112, unit: "%",
              bars: mc.e10Bars,
            },
          ],
        };
      }

      // Vista instantánea de ablación E10 sobre la configuración actual (sin congelar sliders)
      const pBase = Object.assign(V.defaults(), {
        visibility: p.visibility, days_dirty: Math.max(p.days_dirty, 0), T_amb: p.T_amb, I_bg: p.I_bg, Rb: p.Rb,
        P_comm_extra_per_led: p.P_comm_extra_per_led,
      });
      const pFilt = Object.assign({}, pBase, { use_blue_filter: true, days_dirty: Math.min(pBase.days_dirty, 7) });
      const pFec = Object.assign({}, pFilt, { use_fec: true });
      const pAdr = Object.assign(V.optimizedGrid(pFec), { adr_enabled: true });

      const rB = ctx.evalConfigInstant(pBase, idxS);
      const rF = ctx.evalConfigInstant(pFilt, idxS);
      const rC = ctx.evalConfigInstant(pFec, idxS);
      const rA = ctx.evalConfigInstant(pAdr, idxS);

      // Wi-Fi + FEC en los mismos índices
      const wFec = Object.assign({}, w, { use_fec: true });
      const rRng = V.rng(17 + ctx.seed);
      let wOk = 0, wGp = 0, hOk = 0, hGp = 0;
      for (let k = 0; k < idxS.length; k++) {
        const i = idxS[k];
        const o = V.wifiPoint(wFec, grid.pts[3 * i], grid.pts[3 * i + 1], 1, ctx.heavy, ctx.workersB, ctx.sfNorm[i], rRng, 25);
        const wGood = o.avail >= 0.99;
        const vGood = rA.ber[k] <= 1e-6;
        if (wGood) wOk++;
        wGp += o.gp;
        if (wGood || vGood) hOk++;
        hGp += Math.max(o.gp, rA.gp[k]);
      }
      const n = idxS.length;
      const stB = statsVLC(rB), stF = statsVLC(rF), stC = statsVLC(rC), stA = statsVLC(rA);
      const covW = 100 * wOk / n, gpW = wGp / n, covH = 100 * hOk / n, gpH = hGp / n;
      const pwFec = V.wifiPcomm(wFec);
      const barsCov = [
        { label: "1. Base\n3×3", value: stB.cov, color: "fail" },
        { label: "2. +Filtro\n+≤7 d", value: stF.cov, color: "hivis" },
        { label: "3. +RS\nFEC", value: stC.cov, color: "safety" },
        { label: "4. +ADR\n4×4", value: stA.cov, color: "ok" },
        { label: "5. Wi-Fi\n+FEC", value: covW, color: "rf" },
        { label: "6. Híbrido\nVLC+RF", value: covH, color: "blue" },
      ];
      const barsGp = [
        { label: "1. Base\n3×3", value: stB.gp, sub: `${fmt(stB.gp / V.D.Pcomm(pBase), 2)} M/J`, color: "fail" },
        { label: "2. +Filtro\n+≤7 d", value: stF.gp, sub: `${fmt(stF.gp / V.D.Pcomm(pFilt), 2)} M/J`, color: "hivis" },
        { label: "3. +RS\nFEC", value: stC.gp, sub: `${fmt(stC.gp / V.D.Pcomm(pFec), 2)} M/J`, color: "safety" },
        { label: "4. +ADR\n4×4", value: stA.gp, sub: `${fmt(stA.gp / V.D.Pcomm(pAdr), 2)} M/J`, color: "ok" },
        { label: "5. Wi-Fi\n+FEC", value: gpW, sub: `${fmt(gpW / pwFec, 2)} M/J`, color: "rf" },
        { label: "6. Híbrido\nVLC+RF", value: gpH, sub: `${fmt(gpH / (V.D.Pcomm(pAdr) + pwFec), 2)} M/J`, color: "blue" },
      ];
      const maxGp = Math.max(p.Rb / 1e6 * 1.2, 2);
      return {
        tabId,
        mcActive: false,
        panels: [
          {
            kind: "bars",
            title: "E10 · Ablación instantánea: Confiabilidad (escena actual)",
            yLabel: "Área confiable [%]", yMin: 0, yMax: 112, unit: "%",
            bars: barsCov,
          },
          {
            kind: "bars",
            title: "E10 · Ablación instantánea: Goodput [Mb/s] y Mbit/J",
            yLabel: "Goodput medio [Mb/s]", yMin: 0, yMax: maxGp, unit: "",
            bars: barsGp,
          },
        ],
      };
    }

    return { tabId, panels: [] };
  }

  // =========================================================================
  // 2. MOTOR DE RENDERIZADO CANVAS 2D (nativo, Retina, claro/oscuro, tooltips)
  // =========================================================================
  function readPalette() {
    const cs = getComputedStyle(document.documentElement);
    const get = (k, fb) => cs.getPropertyValue(k).trim() || fb;
    const isDark = cs.colorScheme.includes("dark") || document.documentElement.dataset.theme === "dark";
    return {
      bg: get("--panel", "#f2f4f3"),
      panel2: get("--panel-2", "#e9eceb"),
      ink: get("--ink", "#1e2a2f"),
      muted: get("--muted", "#56646a"),
      line: get("--line", "#c3cac9"),
      ok: get("--ok", "#0f7a74"),
      rf: get("--rf", "#b4236d"),
      safety: get("--safety", "#d99b00"),
      hivis: get("--hivis", "#e8641b"),
      fail: get("--fail", "#c0362c"),
      blue: isDark ? "#6aa8ff" : "#1f6fd1",
      purple: isDark ? "#b794f4" : "#6b46c1",
      isDark,
    };
  }

  function setupCanvas(canvas) {
    const rect = canvas.getBoundingClientRect();
    const w = Math.max(220, Math.floor(rect.width));
    const h = Math.max(150, Math.floor(rect.height));
    const dpr = Math.min(window.devicePixelRatio || 1, 2);
    if (canvas.width !== Math.floor(w * dpr) || canvas.height !== Math.floor(h * dpr)) {
      canvas.width = Math.floor(w * dpr);
      canvas.height = Math.floor(h * dpr);
    }
    const ctx = canvas.getContext("2d");
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    return { ctx, w, h };
  }

  function niceTicks(min, max, maxTicks = 5) {
    const span = max - min;
    if (!(span > 0)) return [min];
    const rough = span / maxTicks;
    const pow = Math.pow(10, Math.floor(Math.log10(rough)));
    const norm = rough / pow;
    const step = (norm <= 1 ? 1 : norm <= 2 ? 2 : norm <= 2.5 ? 2.5 : norm <= 5 ? 5 : 10) * pow;
    const start = Math.ceil(min / step - 1e-9) * step;
    const out = [];
    for (let v = start; v <= max + 1e-9; v += step) out.push(+v.toFixed(6));
    return out;
  }

  function drawXYChart(canvas, spec, pal, hoverX) {
    const { ctx, w, h } = setupCanvas(canvas);
    ctx.clearRect(0, 0, w, h);

    const pad = { l: 48, r: 16, t: 38, b: 36 };
    const pw = Math.max(60, w - pad.l - pad.r);
    const ph = Math.max(60, h - pad.t - pad.b);

    // Escalas X
    let xMap, xInv, xTicks;
    if (spec.xScale === "log-inv") {
      const lMax = Math.log10(200), lMin = Math.log10(5);
      xMap = (x) => pad.l + ((lMax - Math.log10(Math.max(x, 1e-6))) / (lMax - lMin)) * pw;
      xInv = (px) => Math.pow(10, lMax - ((px - pad.l) / pw) * (lMax - lMin));
      xTicks = spec.xTicks || [200, 100, 50, 20, 10, 5];
    } else if (spec.xScale === "log") {
      const lMin = Math.log10(1), lMax = Math.log10(100);
      xMap = (x) => pad.l + ((Math.log10(Math.max(x, 1e-6)) - lMin) / (lMax - lMin)) * pw;
      xInv = (px) => Math.pow(10, lMin + ((px - pad.l) / pw) * (lMax - lMin));
      xTicks = spec.xTicks || [1, 2, 5, 10, 20, 50, 100];
    } else {
      const xMin = spec.xMin ?? 0, xMax = spec.xMax ?? 100;
      xMap = (x) => pad.l + ((x - xMin) / Math.max(xMax - xMin, 1e-9)) * pw;
      xInv = (px) => xMin + ((px - pad.l) / pw) * (xMax - xMin);
      xTicks = spec.xTicks || niceTicks(xMin, xMax, 6);
    }

    const yMin = spec.yMin ?? 0, yMax = spec.yMax ?? 100;
    const yMap = (y) => pad.t + ph - ((Math.min(Math.max(y, yMin), yMax) - yMin) / Math.max(yMax - yMin, 1e-9)) * ph;
    const yTicks = spec.yTicks || ((yMin === 0 && yMax >= 100 && yMax <= 115) ? [0, 25, 50, 75, 100]
      : (yMin === -30 && yMax === 0) ? [-30, -24, -18, -12, -6, 0] : niceTicks(yMin, yMax, 5));

    // Fondo del área de trazado
    ctx.fillStyle = pal.panel2;
    ctx.fillRect(pad.l, pad.t, pw, ph);

    // Grilla y ticks
    ctx.strokeStyle = pal.line;
    ctx.lineWidth = 1;
    ctx.font = '11px "Barlow", sans-serif';
    ctx.fillStyle = pal.muted;

    // Horizontales
    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (const yv of yTicks) {
      const py = yMap(yv);
      ctx.globalAlpha = 0.45;
      ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(pad.l + pw, py); ctx.stroke();
      ctx.globalAlpha = 1;
      const dec = Math.abs(yv - Math.round(yv)) > 1e-4 ? 1 : 0;
      ctx.fillText(fmt(yv, dec), pad.l - 6, py);
    }

    // Verticales
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    for (const xv of xTicks) {
      const px = xMap(xv);
      if (px < pad.l - 2 || px > pad.l + pw + 2) continue;
      ctx.globalAlpha = 0.45;
      ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke();
      ctx.globalAlpha = 1;
      const dec = Math.abs(xv - Math.round(xv)) > 1e-4 ? 1 : 0;
      ctx.fillText(fmt(xv, dec), px, pad.t + ph + 5);
    }

    // Borde del área
    ctx.strokeStyle = pal.line; ctx.lineWidth = 1;
    ctx.strokeRect(pad.l, pad.t, pw, ph);

    // Líneas de referencia horizontales y verticales
    if (spec.hLines) {
      for (const hl of spec.hLines) {
        if (hl.y < yMin || hl.y > yMax) continue;
        const py = yMap(hl.y);
        ctx.save();
        ctx.strokeStyle = pal[hl.color] || pal.muted;
        ctx.lineWidth = hl.dashed ? 1.2 : 1.7;
        if (hl.dashed) ctx.setLineDash([4, 4]);
        ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(pad.l + pw, py); ctx.stroke();
        ctx.fillStyle = pal[hl.color] || pal.muted;
        ctx.font = '600 10.5px "Barlow", sans-serif';
        ctx.textAlign = "right"; ctx.textBaseline = "bottom";
        if (hl.label) ctx.fillText(hl.label, pad.l + pw - 4, py - 2);
        ctx.restore();
      }
    }
    if (spec.vLines) {
      for (const vl of spec.vLines) {
        const px = xMap(vl.x);
        if (px < pad.l || px > pad.l + pw) continue;
        ctx.save();
        ctx.strokeStyle = pal[vl.color] || pal.muted;
        ctx.lineWidth = 1.3;
        if (vl.dashed) ctx.setLineDash([4, 4]);
        ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke();
        ctx.fillStyle = pal[vl.color] || pal.muted;
        ctx.font = '600 10.5px "Barlow", sans-serif';
        const leftSide = px > pad.l + pw * 0.45;
        ctx.textAlign = leftSide ? "right" : "left";
        ctx.textBaseline = "bottom";
        if (vl.label) ctx.fillText(vl.label, px + (leftSide ? -4 : 4), pad.t + ph - 4);
        ctx.restore();
      }
    }

    // Anotaciones discretas (ej. regímenes EMI)
    if (spec.annotations) {
      ctx.save();
      ctx.fillStyle = pal.muted;
      ctx.font = '10px "Barlow", sans-serif';
      ctx.textAlign = "center"; ctx.textBaseline = "bottom";
      for (const a of spec.annotations) ctx.fillText(a.text, xMap(a.x), yMap(a.y));
      ctx.restore();
    }

    // Series (clip dentro del área de trazado)
    ctx.save();
    ctx.beginPath(); ctx.rect(pad.l, pad.t, pw, ph); ctx.clip();
    for (const s of spec.series) {
      const col = pal[s.color] || pal.ok;
      ctx.strokeStyle = col;
      ctx.fillStyle = col;
      ctx.lineWidth = s.width || 2;
      ctx.setLineDash(s.dashed ? [5, 4] : []);
      ctx.beginPath();
      let started = false;
      for (let i = 0; i < s.xs.length; i++) {
        if (!Number.isFinite(s.ys[i])) continue;
        const px = xMap(s.xs[i]), py = yMap(s.ys[i]);
        if (!started) { ctx.moveTo(px, py); started = true; } else ctx.lineTo(px, py);
      }
      ctx.stroke();
      ctx.setLineDash([]);
      if (s.dots) {
        for (let i = 0; i < s.xs.length; i++) {
          if (!Number.isFinite(s.ys[i])) continue;
          ctx.beginPath(); ctx.arc(xMap(s.xs[i]), yMap(s.ys[i]), 2.7, 0, Math.PI * 2); ctx.fill();
        }
      }
    }
    ctx.restore();

    // Marcador del valor actual del slider
    if (spec.currentMarker) {
      const cm = spec.currentMarker;
      ctx.save();
      if (Number.isFinite(cm.x)) {
        const px = xMap(cm.x);
        if (px >= pad.l && px <= pad.l + pw) {
          ctx.strokeStyle = pal.safety;
          ctx.lineWidth = 1.6;
          ctx.setLineDash([3, 3]);
          ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke();
          ctx.setLineDash([]);
          // Pastilla superior
          ctx.font = '600 10.5px "Barlow", sans-serif';
          const tw = ctx.measureText(cm.label).width + 10;
          const bx = Math.min(Math.max(px - tw / 2, pad.l + 2), pad.l + pw - tw - 2);
          ctx.fillStyle = pal.safety;
          ctx.fillRect(bx, pad.t + 3, tw, 16);
          ctx.fillStyle = "#1e2a2f";
          ctx.textAlign = "center"; ctx.textBaseline = "middle";
          ctx.fillText(cm.label, bx + tw / 2, pad.t + 11);
        }
      } else if (cm.note) {
        ctx.font = '600 10.5px "Barlow", sans-serif';
        const tw = ctx.measureText(cm.note).width + 10;
        ctx.fillStyle = pal.panel;
        ctx.fillRect(pad.l + 4, pad.t + 4, tw, 16);
        ctx.fillStyle = pal.muted;
        ctx.textAlign = "left"; ctx.textBaseline = "middle";
        ctx.fillText(cm.note, pad.l + 9, pad.t + 12);
      }
      ctx.restore();
    }

    // Guía vertical de hover
    if (hoverX != null && hoverX >= pad.l && hoverX <= pad.l + pw) {
      ctx.save();
      ctx.strokeStyle = pal.ink; ctx.globalAlpha = 0.45; ctx.lineWidth = 1;
      ctx.beginPath(); ctx.moveTo(hoverX, pad.t); ctx.lineTo(hoverX, pad.t + ph); ctx.stroke();
      ctx.restore();
    }

    // Título y leyenda compacta arriba
    ctx.fillStyle = pal.ink;
    ctx.font = '600 13px "Barlow Condensed", sans-serif';
    ctx.textAlign = "left"; ctx.textBaseline = "top";
    ctx.fillText(spec.title, pad.l, 4);

    let lx = pad.l;
    ctx.font = '10.5px "Barlow", sans-serif';
    ctx.textBaseline = "middle";
    for (const s of spec.series) {
      const col = pal[s.color] || pal.ok;
      const lw = ctx.measureText(s.label).width;
      if (lx + lw + 24 > w - 4) break;
      ctx.strokeStyle = col; ctx.lineWidth = 2.2;
      ctx.setLineDash(s.dashed ? [4, 3] : []);
      ctx.beginPath(); ctx.moveTo(lx, 26); ctx.lineTo(lx + 12, 26); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = pal.muted;
      ctx.fillText(s.label, lx + 16, 26);
      lx += lw + 26;
    }

    // Etiquetas de ejes
    ctx.fillStyle = pal.muted;
    ctx.font = '11px "Barlow", sans-serif';
    ctx.textAlign = "center"; ctx.textBaseline = "bottom";
    ctx.fillText(spec.xLabel, pad.l + pw / 2, h - 2);
    ctx.save();
    ctx.translate(12, pad.t + ph / 2);
    ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    ctx.fillText(spec.yLabel, 0, 0);
    ctx.restore();

    canvas._meta = { kind: "xy", pad, pw, ph, xMap, xInv, yMap, spec };
  }

  function drawImpulseChart(canvas, spec, pal, hoverX) {
    const { ctx, w, h } = setupCanvas(canvas);
    ctx.clearRect(0, 0, w, h);
    const pad = { l: 52, r: 16, t: 38, b: 36 };
    const pw = Math.max(60, w - pad.l - pad.r);
    const ph = Math.max(60, h - pad.t - pad.b);
    const xMin = spec.xMin, xMax = spec.xMax;
    const yMax = spec.maxH * 1.22;
    const xMap = (x) => pad.l + ((x - xMin) / (xMax - xMin)) * pw;
    const xInv = (px) => xMin + ((px - pad.l) / pw) * (xMax - xMin);
    const yMap = (y) => pad.t + ph - (Math.min(Math.max(y, 0), yMax) / yMax) * ph;

    ctx.fillStyle = pal.panel2; ctx.fillRect(pad.l, pad.t, pw, ph);
    ctx.strokeStyle = pal.line; ctx.lineWidth = 1;
    ctx.font = '11px "Barlow", sans-serif'; ctx.fillStyle = pal.muted;

    ctx.textAlign = "right"; ctx.textBaseline = "middle";
    for (let i = 0; i <= 4; i++) {
      const yv = yMax * i / 4, py = yMap(yv);
      ctx.globalAlpha = 0.45; ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(pad.l + pw, py); ctx.stroke();
      ctx.globalAlpha = 1; ctx.fillText(yv >= 10 ? fmt(yv, 0) : fmt(yv, 1), pad.l - 6, py);
    }
    ctx.textAlign = "center"; ctx.textBaseline = "top";
    for (let t = 10; t <= 90; t += 10) {
      const px = xMap(t);
      ctx.globalAlpha = 0.45; ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke();
      ctx.globalAlpha = 1; ctx.fillText(String(t), px, pad.t + ph + 5);
    }
    ctx.strokeStyle = pal.line; ctx.strokeRect(pad.l, pad.t, pw, ph);

    // Curva de reflexiones con relleno suave
    ctx.save();
    ctx.beginPath(); ctx.rect(pad.l, pad.t, pw, ph); ctx.clip();
    if (spec.tNs.length) {
      ctx.beginPath();
      ctx.moveTo(xMap(spec.tNs[0]), pad.t + ph);
      for (let i = 0; i < spec.tNs.length; i++) ctx.lineTo(xMap(spec.tNs[i]), yMap(spec.hVals[i]));
      ctx.lineTo(xMap(spec.tNs[spec.tNs.length - 1]), pad.t + ph);
      ctx.closePath();
      ctx.fillStyle = pal.blue; ctx.globalAlpha = 0.16; ctx.fill();
      ctx.globalAlpha = 1;
      ctx.beginPath();
      for (let i = 0; i < spec.tNs.length; i++) {
        const px = xMap(spec.tNs[i]), py = yMap(spec.hVals[i]);
        if (i === 0) ctx.moveTo(px, py); else ctx.lineTo(px, py);
      }
      ctx.strokeStyle = pal.blue; ctx.lineWidth = 2; ctx.stroke();
    }

    // Llegadas LOS verticales
    let maxG = 0;
    for (const a of spec.losArrivals) if (a.g > maxG) maxG = a.g;
    for (const a of spec.losArrivals) {
      const px = xMap(a.tNs);
      const frac = maxG > 0 ? 0.35 + 0.60 * (a.g / maxG) : 0.8;
      const pyTop = pad.t + ph * (1 - frac);
      ctx.strokeStyle = pal.ok; ctx.lineWidth = 2.4;
      ctx.beginPath(); ctx.moveTo(px, pad.t + ph); ctx.lineTo(px, pyTop); ctx.stroke();
      ctx.fillStyle = pal.ok;
      ctx.beginPath(); ctx.arc(px, pyTop, 3.2, 0, Math.PI * 2); ctx.fill();
    }

    // Instante óptimo de muestreo t0
    if (Number.isFinite(spec.t0Ns)) {
      const px = xMap(spec.t0Ns);
      ctx.strokeStyle = pal.safety; ctx.lineWidth = 1.6; ctx.setLineDash([4, 3]);
      ctx.beginPath(); ctx.moveTo(px, pad.t); ctx.lineTo(px, pad.t + ph); ctx.stroke();
      ctx.setLineDash([]);
      ctx.fillStyle = pal.safety;
      const lbl = `t₀ ópt. ${fmt(spec.t0Ns, 1)} ns`;
      ctx.font = '600 10.5px "Barlow", sans-serif';
      const tw = ctx.measureText(lbl).width + 8;
      const bx = Math.min(Math.max(px - tw / 2, pad.l + 2), pad.l + pw - tw - 2);
      ctx.fillRect(bx, pad.t + 3, tw, 16);
      ctx.fillStyle = "#1e2a2f"; ctx.textAlign = "center"; ctx.textBaseline = "middle";
      ctx.fillText(lbl, bx + tw / 2, pad.t + 11);
    }
    ctx.restore();

    if (hoverX != null && hoverX >= pad.l && hoverX <= pad.l + pw) {
      ctx.save(); ctx.strokeStyle = pal.ink; ctx.globalAlpha = 0.45;
      ctx.beginPath(); ctx.moveTo(hoverX, pad.t); ctx.lineTo(hoverX, pad.t + ph); ctx.stroke();
      ctx.restore();
    }

    ctx.fillStyle = pal.ink; ctx.font = '600 13px "Barlow Condensed", sans-serif';
    ctx.textAlign = "left"; ctx.textBaseline = "top";
    ctx.fillText(spec.title, pad.l, 4);

    // Leyenda compacta
    ctx.font = '10.5px "Barlow", sans-serif'; ctx.textBaseline = "middle";
    ctx.fillStyle = pal.ok; ctx.fillRect(pad.l, 24, 10, 4);
    ctx.fillStyle = pal.muted; ctx.fillText("Llegadas LOS (por luminaria)", pad.l + 14, 26);
    ctx.fillStyle = pal.blue; ctx.fillRect(pad.l + 165, 24, 10, 4);
    ctx.fillStyle = pal.muted; ctx.fillText("Reflexiones en paredes h(t)", pad.l + 179, 26);

    ctx.fillStyle = pal.muted; ctx.font = '11px "Barlow", sans-serif';
    ctx.textAlign = "center"; ctx.textBaseline = "bottom"; ctx.fillText(spec.xLabel, pad.l + pw / 2, h - 2);
    ctx.save(); ctx.translate(12, pad.t + ph / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillText(spec.yLabel, 0, 0); ctx.restore();

    canvas._meta = { kind: "impulse", pad, pw, ph, xMap, xInv, spec };
  }

  function drawBarChart(canvas, spec, pal, hoverX) {
    const { ctx, w, h } = setupCanvas(canvas);
    ctx.clearRect(0, 0, w, h);
    const isGroupedLog = spec.kind === "bars-grouped-log";
    const pad = { l: isGroupedLog ? 52 : 44, r: 12, t: 38, b: 36 };
    const pw = Math.max(60, w - pad.l - pad.r);
    const ph = Math.max(60, h - pad.t - pad.b);

    ctx.fillStyle = pal.panel2; ctx.fillRect(pad.l, pad.t, pw, ph);
    ctx.strokeStyle = pal.line; ctx.lineWidth = 1;
    ctx.font = '11px "Barlow", sans-serif'; ctx.fillStyle = pal.muted;

    if (isGroupedLog) {
      const eMin = spec.yMinExp, eMax = spec.yMaxExp;
      const yMap = (val) => {
        const lv = Math.min(Math.max(Math.log10(Math.max(val, Math.pow(10, eMin))), eMin), eMax);
        return pad.t + ph - ((lv - eMin) / (eMax - eMin)) * ph;
      };
      ctx.textAlign = "right"; ctx.textBaseline = "middle";
      for (let e = eMin; e <= eMax; e++) {
        if ((eMax - eMin) > 5 && (e - eMin) % 2 !== 0 && e !== eMax) continue;
        const py = pad.t + ph - ((e - eMin) / (eMax - eMin)) * ph;
        ctx.globalAlpha = 0.45; ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(pad.l + pw, py); ctx.stroke();
        ctx.globalAlpha = 1; ctx.fillText(`10^${e}`, pad.l - 5, py);
      }
      ctx.strokeRect(pad.l, pad.t, pw, ph);

      const nG = spec.labels.length, nS = spec.groups.length;
      const slot = pw / nG, bw = Math.min(22, (slot * 0.72) / nS);
      for (let g = 0; g < nG; g++) {
        const cx = pad.l + (g + 0.5) * slot;
        if (g === spec.activeGroup) {
          ctx.fillStyle = pal.safety; ctx.globalAlpha = 0.12;
          ctx.fillRect(pad.l + g * slot + 2, pad.t + 1, slot - 4, ph - 2);
          ctx.globalAlpha = 1;
        }
        const startX = cx - (nS * bw) / 2;
        for (let s = 0; s < nS; s++) {
          const val = spec.groups[s].values[g];
          const py = yMap(val);
          ctx.fillStyle = pal[spec.groups[s].color] || pal.ok;
          ctx.fillRect(startX + s * bw + 1, py, bw - 2, pad.t + ph - py);
        }
        ctx.fillStyle = g === spec.activeGroup ? pal.ink : pal.muted;
        ctx.font = (g === spec.activeGroup ? "600 " : "") + '10.5px "Barlow", sans-serif';
        ctx.textAlign = "center"; ctx.textBaseline = "top";
        spec.labels[g].split("\n").forEach((ln, li) => ctx.fillText(ln, cx, pad.t + ph + 4 + li * 12));
      }

      // Leyenda
      let lx = pad.l; ctx.font = '10.5px "Barlow", sans-serif'; ctx.textBaseline = "middle"; ctx.textAlign = "left";
      for (const s of spec.groups) {
        const lw = ctx.measureText(s.label).width;
        ctx.fillStyle = pal[s.color] || pal.ok; ctx.fillRect(lx, 23, 9, 7);
        ctx.fillStyle = pal.muted; ctx.fillText(s.label, lx + 12, 26);
        lx += lw + 22;
      }
      canvas._meta = { kind: "bars-grouped-log", pad, pw, ph, slot, spec };
    } else {
      const yMin = spec.yMin ?? 0, yMax = spec.yMax ?? 100;
      const yMap = (y) => pad.t + ph - ((Math.min(Math.max(y, yMin), yMax) - yMin) / Math.max(yMax - yMin, 1e-9)) * ph;
      const bTicks = (yMin === 0 && yMax >= 100 && yMax <= 115) ? [0, 25, 50, 75, 100] : niceTicks(yMin, yMax, 4);
      ctx.textAlign = "right"; ctx.textBaseline = "middle";
      for (const yv of bTicks) {
        const py = yMap(yv);
        ctx.globalAlpha = 0.45; ctx.beginPath(); ctx.moveTo(pad.l, py); ctx.lineTo(pad.l + pw, py); ctx.stroke();
        ctx.globalAlpha = 1; ctx.fillText(fmt(yv, Math.abs(yv - Math.round(yv)) > 1e-4 ? 1 : 0), pad.l - 5, py);
      }
      ctx.strokeRect(pad.l, pad.t, pw, ph);

      const nB = spec.bars.length, slot = pw / nB, bw = Math.min(44, slot * 0.64);
      const hoverIdx = (hoverX != null && hoverX >= pad.l && hoverX <= pad.l + pw) ? Math.floor((hoverX - pad.l) / slot) : -1;
      for (let i = 0; i < nB; i++) {
        const b = spec.bars[i], cx = pad.l + (i + 0.5) * slot, py = yMap(b.value);
        ctx.fillStyle = pal[b.color] || pal.ok;
        ctx.globalAlpha = hoverIdx === i ? 1 : 0.88;
        ctx.fillRect(cx - bw / 2, py, bw, pad.t + ph - py);
        ctx.globalAlpha = 1;
        ctx.fillStyle = pal.ink;
        ctx.font = '600 10.5px "Barlow", sans-serif'; ctx.textAlign = "center"; ctx.textBaseline = "bottom";
        const valTxt = spec.unit === "%" ? `${fmt(b.value, 1)}%` : fmt(b.value, 2);
        ctx.fillText(valTxt, cx, py - (b.sub ? 12 : 3));
        if (b.sub) {
          ctx.fillStyle = pal.muted; ctx.font = '9.5px "Barlow", sans-serif';
          ctx.fillText(b.sub, cx, py - 2);
        }
        ctx.fillStyle = pal.muted; ctx.font = '10px "Barlow", sans-serif'; ctx.textBaseline = "top";
        b.label.split("\n").forEach((ln, li) => ctx.fillText(ln, cx, pad.t + ph + 4 + li * 12));
      }
      canvas._meta = { kind: "bars", pad, pw, ph, slot, spec };
    }

    ctx.fillStyle = pal.ink; ctx.font = '600 13px "Barlow Condensed", sans-serif';
    ctx.textAlign = "left"; ctx.textBaseline = "top";
    ctx.fillText(spec.title, pad.l, 4);

    ctx.fillStyle = pal.muted; ctx.font = '11px "Barlow", sans-serif';
    ctx.save(); ctx.translate(11, pad.t + ph / 2); ctx.rotate(-Math.PI / 2);
    ctx.textAlign = "center"; ctx.textBaseline = "top"; ctx.fillText(spec.yLabel, 0, 0); ctx.restore();
  }

  function tooltipHtml(meta, localX) {
    if (!meta || localX < meta.pad.l || localX > meta.pad.l + meta.pw) return null;
    if (meta.kind === "xy") {
      const xv = meta.xInv(localX);
      const spec = meta.spec;
      const lines = [`<b>${spec.xLabel.split("[")[0].trim()}: ${fmt(xv, xv < 10 ? 2 : 1)}</b>`];
      for (const s of spec.series) {
        if (!s.xs.length) continue;
        let best = 0, bd = Infinity;
        for (let i = 0; i < s.xs.length; i++) {
          const d = Math.abs(meta.xMap(s.xs[i]) - localX);
          if (d < bd) { bd = d; best = i; }
        }
        lines.push(`${s.label}: <b>${fmt(s.ys[best], 2)}</b>`);
      }
      return lines.join("<br>");
    }
    if (meta.kind === "impulse") {
      const t = meta.xInv(localX);
      const spec = meta.spec;
      let best = 0, bd = Infinity;
      for (let i = 0; i < spec.tNs.length; i++) {
        const d = Math.abs(spec.tNs[i] - t);
        if (d < bd) { bd = d; best = i; }
      }
      const nearLos = spec.losArrivals.filter((a) => Math.abs(a.tNs - t) < 1.2);
      const lines = [`<b>t = ${fmt(t, 1)} ns</b>`, `h_refl(t): <b>${fmt(spec.hVals[best] || 0, 1)} 1/s</b>`];
      if (nearLos.length) lines.push(`LOS lum. #${nearLos.map((a) => a.j + 1).join(", #")}`);
      return lines.join("<br>");
    }
    if (meta.kind === "bars-grouped-log") {
      const g = Math.min(meta.spec.labels.length - 1, Math.max(0, Math.floor((localX - meta.pad.l) / meta.slot)));
      const lines = [`<b>${meta.spec.labels[g].replace("\n", " · ")}</b>`];
      for (const s of meta.spec.groups) lines.push(`${s.label}: <b>${sciText(s.values[g])} A²</b>`);
      return lines.join("<br>");
    }
    if (meta.kind === "bars") {
      const idx = Math.min(meta.spec.bars.length - 1, Math.max(0, Math.floor((localX - meta.pad.l) / meta.slot)));
      const b = meta.spec.bars[idx];
      const vTxt = meta.spec.unit === "%" ? `${fmt(b.value, 1)} %` : `${fmt(b.value, 2)} Mb/s`;
      return `<b>${b.label.replace("\n", " ")}</b><br>${vTxt}${b.sub ? ` (${b.sub})` : ""}`;
    }
    return null;
  }

  function renderPanels(hostEl, tipEl, tabData) {
    const pal = readPalette();
    const panels = tabData.panels || [];
    hostEl.dataset.cols = String(panels.length > 1 ? 2 : 1);

    // Reutilizar o crear los <canvas>
    while (hostEl.children.length > panels.length) hostEl.removeChild(hostEl.lastChild);
    while (hostEl.children.length < panels.length) {
      const cv = document.createElement("canvas");
      cv.className = "chart-cv";
      cv.addEventListener("pointermove", (ev) => {
        const r = cv.getBoundingClientRect();
        const lx = ev.clientX - r.left;
        const pSpec = cv._panelSpec;
        if (!pSpec) return;
        const palNow = readPalette();
        if (pSpec.kind === "xy") drawXYChart(cv, pSpec, palNow, lx);
        else if (pSpec.kind === "impulse") drawImpulseChart(cv, pSpec, palNow, lx);
        else drawBarChart(cv, pSpec, palNow, lx);
        const html = tooltipHtml(cv._meta, lx);
        if (html) {
          tipEl.innerHTML = html;
          tipEl.hidden = false;
          const hr = hostEl.getBoundingClientRect();
          const tx = Math.min(Math.max(ev.clientX - hr.left + 12, 6), hr.width - tipEl.offsetWidth - 8);
          const ty = Math.min(Math.max(ev.clientY - hr.top - 10, 6), hr.height - tipEl.offsetHeight - 8);
          tipEl.style.left = tx + "px";
          tipEl.style.top = ty + "px";
        } else tipEl.hidden = true;
      });
      cv.addEventListener("pointerleave", () => {
        tipEl.hidden = true;
        const pSpec = cv._panelSpec;
        if (!pSpec) return;
        const palNow = readPalette();
        if (pSpec.kind === "xy") drawXYChart(cv, pSpec, palNow, null);
        else if (pSpec.kind === "impulse") drawImpulseChart(cv, pSpec, palNow, null);
        else drawBarChart(cv, pSpec, palNow, null);
      });
      hostEl.appendChild(cv);
    }

    panels.forEach((pSpec, i) => {
      const cv = hostEl.children[i];
      cv._panelSpec = pSpec;
      if (pSpec.kind === "xy") drawXYChart(cv, pSpec, pal, null);
      else if (pSpec.kind === "impulse") drawImpulseChart(cv, pSpec, pal, null);
      else drawBarChart(cv, pSpec, pal, null);
    });
  }

  return { computeTabData, renderPanels, sampleIndices, statsVLC };
})();
if (typeof module !== "undefined") module.exports = VLCCharts;
