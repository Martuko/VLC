/* vlc_physics.js — port a JavaScript del núcleo vlc_core.py (modelo de ISI por apertura de ojo).
 * Sin dependencias ni DOM: se usa en el navegador (inline en index.html) y en Node para validar contra Python.
 * Unidades SI. Simplificación declarada: los obstáculos móviles bloquean la línea de vista, no las reflexiones. */
const VLC = (() => {
  "use strict";
  const Q_E = 1.602176634e-19, K_B = 1.380649e-23, C = 2.99792458e8;

  // ---------------------------------------------------------------- funciones especiales
  function erfc(x) { // Numerical Recipes erfcc, error relativo < 1.2e-7 en todo el dominio
    const z = Math.abs(x), t = 1 / (1 + 0.5 * z);
    const r = t * Math.exp(-z * z - 1.26551223 + t * (1.00002368 + t * (0.37409196 + t * (0.09678418 +
      t * (-0.18628806 + t * (0.27886807 + t * (-1.13520398 + t * (1.48851587 + t * (-0.82215223 + t * 0.17087277)))))))));
    return x >= 0 ? r : 2 - r;
  }
  const qfunc = (x) => 0.5 * erfc(x / Math.SQRT2);
  function lgamma(x) { // Lanczos
    const g = [76.18009172947146, -86.50532032941677, 24.01409824083091, -1.231739572450155,
      0.1208650973866179e-2, -0.5395239384953e-5];
    let y = x, tmp = x + 5.5; tmp -= (x + 0.5) * Math.log(tmp);
    let ser = 1.000000000190015;
    for (let j = 0; j < 6; j++) ser += g[j] / ++y;
    return -tmp + Math.log(2.5066282746310005 * ser / x);
  }
  function betacf(a, b, x) {
    const MAXIT = 300, EPS = 3e-16, FPMIN = 1e-300;
    let qab = a + b, qap = a + 1, qam = a - 1, c = 1, d = 1 - qab * x / qap;
    if (Math.abs(d) < FPMIN) d = FPMIN; d = 1 / d; let h = d;
    for (let m = 1; m <= MAXIT; m++) {
      const m2 = 2 * m;
      let aa = m * (b - m) * x / ((qam + m2) * (a + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; h *= d * c;
      aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2));
      d = 1 + aa * d; if (Math.abs(d) < FPMIN) d = FPMIN; c = 1 + aa / c; if (Math.abs(c) < FPMIN) c = FPMIN;
      d = 1 / d; const del = d * c; h *= del;
      if (Math.abs(del - 1) < EPS) break;
    }
    return h;
  }
  function betai(a, b, x) { // función beta incompleta regularizada I_x(a,b)
    if (x <= 0) return 0; if (x >= 1) return 1;
    const bt = Math.exp(lgamma(a + b) - lgamma(a) - lgamma(b) + a * Math.log(x) + b * Math.log1p(-x));
    return x < (a + 1) / (a + b + 2) ? bt * betacf(a, b, x) / a : 1 - bt * betacf(b, a, 1 - x) / b;
  }
  function rsPostFecBer(pb, n = 255, k = 239, m = 8) {
    const t = (n - k) >> 1;
    pb = Math.min(Math.max(pb, 0), 0.5);
    const ps = Math.min(Math.max(-Math.expm1(m * Math.log1p(-pb)), 0), 1);
    return pb * betai(t, n - t, ps);
  }
  function goodputMbps(Rb, ber, rate = 1, L = 1000) {
    const pb = Math.min(Math.max(ber, 0), 1 - 1e-15);
    const per = Math.min(Math.max(-Math.expm1(L * Math.log1p(-pb)), 0), 1);
    return (Rb / 1e6) * rate * (1 - per);
  }
  function hgFovFraction(fov, g) {
    const c = Math.cos(Math.min(Math.max(fov, 0), Math.PI));
    if (Math.abs(g) < 1e-6) return 0.5 * (1 - c);
    const v = ((1 - g * g) / (2 * g)) * (1 / (1 - g) - 1 / Math.sqrt(1 + g * g - 2 * g * c));
    return Math.min(Math.max(v, 0), 1);
  }
  function ledBwEff(Rb, f3db) {
    if (!isFinite(f3db)) return 1;
    const tau = 1 / (2 * Math.PI * f3db);
    return Math.min(Math.max(-Math.expm1(-(1 / Rb) / tau), 0), 1);
  }

  // ---------------------------------------------------------------- parámetros (idénticos a Params de vlc_core.py)
  const GRID3 = [[2.5, 2.5], [2.5, 7.5], [2.5, 12.5], [7.5, 2.5], [7.5, 7.5], [7.5, 12.5], [12.5, 2.5], [12.5, 7.5], [12.5, 12.5]];
  const C4 = [1.1, 5.4, 9.6, 13.9];
  const GRID4 = []; for (const x of C4) for (const y of C4) GRID4.push([x, y]);
  function defaults() {
    return {
      L: 15, W: 15, H: 6, rho_wall: 0.4, dA_wall: 0.75,
      led_xy: GRID3, z_led: 5.8, P_opt: 45, led_wall_plug_eff: 0.45, P_comm_extra_per_led: 1.0, P_rx_vlc: 0.5,
      LER: 300, phi_half_deg: 60, mod_index: 0.2,
      f_3db_led: 15e6, use_blue_filter: false, f_3db_blue: 25e6, Ts_blue_signal: 0.55, Ts_blue_bg: 0.25,
      z_rx: 1.0, A_pd: 1e-4, fov_deg: 60, n_conc: 1.5, Ts_filter: 1.0, resp: 0.53,
      adr_enabled: false, adr_tilt_deg: 35, adr_fov_deg: 45,
      I_bg: 1000e-6, I2: 0.562, I3: 0.0868, G_ol: 10, eta_cap: 112e-12 / 1e-4, Gamma: 1.5, g_m: 30e-3, I_dark0: 10e-9,
      Rb: 10e6, use_fec: false, packet_bits: 1000,
      T_amb: 25, delta_T_j_nom: 30.25, T_j_ref: 55.25, k_T_led: 0.004, T_j_derate: 85, derate_slope: 0.015, k_soil_th: 0.25,
      visibility: Infinity, dust_albedo: 0.65, dust_g_hg: 0.70, enable_scattering: true,
      days_dirty: 0, tau_rx_days: 30, tau_tx_days: 120,
    };
  }
  function optimizedGrid(p) {
    const n0 = p.led_xy.length;
    return Object.assign({}, p, { led_xy: GRID4, P_opt: n0 * p.P_opt / 16, phi_half_deg: 65,
      P_comm_extra_per_led: p.P_comm_extra_per_led * n0 / 16 });
  }
  const D = {
    m: (p) => -Math.log(2) / Math.log(Math.cos(p.phi_half_deg * Math.PI / 180)),
    beta: (p) => (isFinite(p.visibility) ? 3.912 / p.visibility : 0),
    Tsrx: (p) => Math.exp(-p.days_dirty / p.tau_rx_days),
    Tstx: (p) => Math.exp(-p.days_dirty / p.tau_tx_days),
    Tj: (p) => p.T_amb + p.delta_T_j_nom * (1 + p.k_soil_th * (1 - Math.exp(-p.days_dirty / p.tau_tx_days))),
    derating(p) {
      const tj = D.Tj(p); let f = 1 - p.k_T_led * (tj - p.T_j_ref);
      if (p.T_j_derate != null && tj > p.T_j_derate) f *= Math.max(0, 1 - p.derate_slope * (tj - p.T_j_derate));
      return Math.max(0, f);
    },
    f3db: (p) => (p.use_blue_filter ? p.f_3db_blue : p.f_3db_led),
    TsBlue: (p) => (p.use_blue_filter ? p.Ts_blue_signal : 1),
    codeRate: (p) => (p.use_fec ? 239 / 255 : 1),
    Pcomm: (p) => p.led_xy.length * p.P_comm_extra_per_led + p.P_rx_vlc * (p.adr_enabled ? 1.6 : 1),
    IbgEff: (p, A) => p.I_bg * (A / 1e-4) * Math.exp(-p.days_dirty / p.tau_rx_days) * p.Ts_filter * (p.use_blue_filter ? p.Ts_blue_bg : 1),
  };

  // ---------------------------------------------------------------- planta
  const STATIC = [
    { name: "Estantería A", b: [1.0, 7.0, 13.2, 14.2, 5.0] },
    { name: "Estantería B", b: [8.0, 14.0, 13.2, 14.2, 5.0] },
    { name: "Prensa", b: [1.5, 3.5, 1.5, 3.5, 4.0] },
    { name: "CNC 1", b: [5.5, 8.5, 2.5, 4.5, 2.5] },
    { name: "CNC 2", b: [10.5, 13.5, 1.5, 4.0, 2.5] },
    { name: "Horno de inducción", b: [10.0, 13.0, 7.0, 10.0, 3.5] },
    { name: "Cinta + soldadura", b: [1.0, 8.5, 8.0, 8.8, 1.3] },
  ];
  const STATIC_BOXES = STATIC.map((s) => s.b);
  const EMI_SOURCES = [
    { name: "Horno de inducción", s: [11.5, 8.5, 2.0, 38.0] },
    { name: "Puesto soldadura", s: [4.75, 8.4, 1.5, 36.0] },
    { name: "Prensa (VFD)", s: [2.5, 2.5, 2.0, 30.0] },
    { name: "CNC 1 (VFD)", s: [7.0, 3.5, 1.8, 28.0] },
    { name: "CNC 2 (VFD)", s: [12.0, 2.75, 1.8, 28.0] },
  ];

  function rng(seed) { // mulberry32
    let a = seed >>> 0;
    const f = () => { a = (a + 0x6D2B79F5) >>> 0; let t = a; t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61); return ((t ^ (t >>> 14)) >>> 0) / 4294967296; };
    f.normal = () => { let u = 0; while (u === 0) u = f(); return Math.sqrt(-2 * Math.log(u)) * Math.cos(2 * Math.PI * f()); };
    return f;
  }
  function randomDynamic(r, nWorkers, nForklifts, p) {
    const out = [];
    const free = (b) => STATIC_BOXES.every((s) => b[1] < s[0] || b[0] > s[1] || b[3] < s[2] || b[2] > s[3]);
    for (const [n, sx, sy, h] of [[nWorkers, 0.5, 0.5, 1.8], [nForklifts, 1.2, 2.4, 2.3]]) {
      let k = 0, tries = 0;
      while (k < n && tries < 50 * Math.max(n, 1)) {
        tries++;
        const rot = r() < 0.5, a = rot ? sy : sx, c = rot ? sx : sy;
        const cx = a / 2 + r() * (p.L - a), cy = c / 2 + r() * (p.W - c);
        const b = [cx - a / 2, cx + a / 2, cy - c / 2, cy + c / 2, h];
        if (free(b)) { out.push(b); k++; }
      }
    }
    return out;
  }
  function segBlocked(px, py, pz, qx, qy, qz, boxes) {
    const eps = 1e-12;
    const dx = qx - px, dy = qy - py, dz = qz - pz;
    const ix = 1 / (Math.abs(dx) < eps ? (dx >= 0 ? eps : -eps) : dx);
    const iy = 1 / (Math.abs(dy) < eps ? (dy >= 0 ? eps : -eps) : dy);
    const iz = 1 / (Math.abs(dz) < eps ? (dz >= 0 ? eps : -eps) : dz);
    for (let k = 0; k < boxes.length; k++) {
      const b = boxes[k];
      let t1 = (b[0] - px) * ix, t2 = (b[1] - px) * ix;
      let tmin = Math.min(t1, t2), tmax = Math.max(t1, t2);
      t1 = (b[2] - py) * iy; t2 = (b[3] - py) * iy;
      tmin = Math.max(tmin, Math.min(t1, t2)); tmax = Math.min(tmax, Math.max(t1, t2));
      t1 = (0 - pz) * iz; t2 = (b[4] - pz) * iz;
      tmin = Math.max(tmin, Math.min(t1, t2)); tmax = Math.min(tmax, Math.max(t1, t2));
      if (tmax >= Math.max(tmin, 0) && tmin <= 1 && tmax > 1e-6) return true;
    }
    return false;
  }
  const insideAny = (x, y, boxes) => boxes.some((b) => x >= b[0] && x <= b[1] && y >= b[2] && y <= b[3]);
  function receiverGrid(p, step) {
    const xs = [], ys = [];
    for (let v = step / 2; v < p.L - 1e-9; v += step) xs.push(+v.toFixed(6));
    for (let v = step / 2; v < p.W - 1e-9; v += step) ys.push(+v.toFixed(6));
    const pts = new Float64Array(xs.length * ys.length * 3); let i = 0;
    for (const y of ys) for (const x of xs) { pts[i++] = x; pts[i++] = y; pts[i++] = p.z_rx; }
    return { pts, xs, ys, n: xs.length * ys.length };
  }
  function walls(p) {
    const d = p.dA_wall, c = [], n = [];
    const zs = [], xs = [], ys = [];
    for (let v = d / 2; v < p.H - 1e-9; v += d) zs.push(v);
    for (let v = d / 2; v < p.L - 1e-9; v += d) xs.push(v);
    for (let v = d / 2; v < p.W - 1e-9; v += d) ys.push(v);
    for (const x of xs) for (const z of zs) { c.push([x, 0, z]); n.push([0, 1, 0]); c.push([x, p.W, z]); n.push([0, -1, 0]); }
    for (const y of ys) for (const z of zs) { c.push([0, y, z]); n.push([1, 0, 0]); c.push([p.L, y, z]); n.push([-1, 0, 0]); }
    return { c, n, dA: d * d };
  }
  function adrNormals(tiltDeg) {
    const th = tiltDeg * Math.PI / 180, s = Math.sin(th), c = Math.cos(th);
    return [[0, 0, 1], [s, 0, c], [0, s, c], [-s, 0, c], [0, -s, c]];
  }
  // Ramas del receptor: 1 fotodiodo vertical o ADR piramidal (5 ramas, área total conservada)
  function branches(p) {
    if (!p.adr_enabled) return [{ n: [0, 0, 1], fov: p.fov_deg * Math.PI / 180, A: p.A_pd }];
    return adrNormals(p.adr_tilt_deg).map((n) => ({ n, fov: p.adr_fov_deg * Math.PI / 180, A: p.A_pd / 5 }));
  }

  // ---------------------------------------------------------------- canal
  /* LOS para un receptor: devuelve, por rama, ganancias h[led] y distancias d[led]. */
  function losChannel(p, x, y, z, boxes, brs) {
    const m = D.m(p), nl = p.led_xy.length;
    const d = new Float64Array(nl), base = new Float64Array(nl), vx = new Float64Array(nl), vy = new Float64Array(nl), vz = new Float64Array(nl);
    for (let j = 0; j < nl; j++) {
      const lx = p.led_xy[j][0], ly = p.led_xy[j][1], lz = p.z_led;
      const ax = lx - x, ay = ly - y, az = lz - z, dd = Math.hypot(ax, ay, az);
      d[j] = dd; vx[j] = ax / dd; vy[j] = ay / dd; vz[j] = az / dd;
      const cphi = az / dd;
      if (cphi <= 0 || (boxes.length && segBlocked(x, y, z, lx, ly, lz, boxes))) { base[j] = 0; continue; }
      base[j] = (m + 1) / (2 * Math.PI * dd * dd) * Math.pow(cphi, m) * p.Ts_filter;
    }
    return brs.map((br) => {
      const g = p.n_conc * p.n_conc / Math.pow(Math.sin(br.fov), 2), h = new Float64Array(nl);
      for (let j = 0; j < nl; j++) {
        if (base[j] === 0) continue;
        const cpsi = vx[j] * br.n[0] + vy[j] * br.n[1] + vz[j] * br.n[2];
        if (cpsi > 0 && Math.acos(Math.min(cpsi, 1)) <= br.fov) h[j] = base[j] * br.A * g * cpsi;
      }
      return { h, d };
    });
  }
  /* Histograma de reflexiones de 1.er orden (sin polvo), por rama: Float32Array(nRx * nb). */
  function reflections(p, grid, boxes, brs, dt = 1e-9, nb = 300, onProgress) {
    const W = walls(p), m = D.m(p), nl = p.led_xy.length, nw = W.c.length;
    const a1 = new Float64Array(nl * nw), d1 = new Float64Array(nl * nw);
    for (let j = 0; j < nl; j++) {
      const lx = p.led_xy[j][0], ly = p.led_xy[j][1], lz = p.z_led;
      for (let w = 0; w < nw; w++) {
        const [cx, cy, cz] = W.c[w], [nx, ny, nz] = W.n[w];
        const vx = cx - lx, vy = cy - ly, vz = cz - lz, dd = Math.hypot(vx, vy, vz);
        const cphi = -vz / dd, calpha = -(vx * nx + vy * ny + vz * nz) / dd;
        d1[j * nw + w] = dd;
        if (cphi > 0 && calpha > 0 && !(boxes.length && segBlocked(cx, cy, cz, lx, ly, lz, boxes)))
          a1[j * nw + w] = (m + 1) / (2 * Math.PI) * Math.pow(cphi, m) * calpha / (dd * dd);
      }
    }
    const out = brs.map(() => new Float32Array(grid.n * nb));
    const gb = brs.map((br) => p.n_conc * p.n_conc / Math.pow(Math.sin(br.fov), 2));
    for (let i = 0; i < grid.n; i++) {
      const x = grid.pts[3 * i], y = grid.pts[3 * i + 1], z = grid.pts[3 * i + 2];
      for (let w = 0; w < nw; w++) {
        const [cx, cy, cz] = W.c[w], [nx, ny, nz] = W.n[w];
        const vx = x - cx, vy = y - cy, vz = z - cz, d2 = Math.hypot(vx, vy, vz);
        const cbeta = (vx * nx + vy * ny + vz * nz) / d2;
        if (cbeta <= 0) continue;
        let blockedChecked = false, blocked = false;
        for (let b = 0; b < brs.length; b++) {
          const br = brs[b];
          const cpsi = -(vx * br.n[0] + vy * br.n[1] + vz * br.n[2]) / d2;
          if (cpsi <= 0 || Math.acos(Math.min(cpsi, 1)) > br.fov) continue;
          if (!blockedChecked) { blocked = boxes.length > 0 && segBlocked(x, y, z, cx, cy, cz, boxes); blockedChecked = true; }
          if (blocked) break;
          const a2 = p.rho_wall * W.dA * cbeta / (Math.PI * d2 * d2) * br.A * p.Ts_filter * gb[b] * cpsi;
          const o = out[b], off = i * nb;
          for (let j = 0; j < nl; j++) {
            const a = a1[j * nw + w]; if (a === 0) continue;
            let k = Math.floor((d1[j * nw + w] + d2) / C / dt); if (k > nb - 1) k = nb - 1;
            o[off + k] += a * a2;
          }
        }
      }
      if (onProgress && (i & 63) === 0) onProgress(i / grid.n);
    }
    return { hist: out, dt, nb };
  }

  // ---------------------------------------------------------------- rendimiento (modelo de ojo)
  function noiseVar(p, Pdc, A) {
    const B = p.Rb, TK = p.T_amb + 273.15;
    const Idark = p.I_dark0 * (A / 1e-4) * Math.pow(2, (p.T_amb - 25) / 10);
    const sShot = 2 * Q_E * p.resp * Pdc * B + 2 * Q_E * D.IbgEff(p, A) * p.I2 * B + 2 * Q_E * Idark * B;
    const Cc = p.eta_cap * A;
    const sTh = (8 * Math.PI * K_B * TK / p.G_ol) * Cc * p.I2 * B * B + (16 * Math.PI * Math.PI * K_B * TK * p.Gamma / p.g_m) * Cc * Cc * p.I3 * B * B * B;
    return { sShot, sTh, v: sShot + sTh };
  }
  const _t = new Float64Array(64), _g = new Float64Array(64), _ord = new Int32Array(64);
  /* Términos de amplitud de una rama: señal, ISI tardía, ISI temprana y varianza de ruido. */
  function eyeTerms(p, br, h, d, hist, off, dt, nb) {
    const beta = D.beta(p), scale = p.P_opt * D.derating(p) * D.Tsrx(p) * D.Tstx(p) * D.TsBlue(p);
    const nl = h.length, Tb = 1 / p.Rb;
    let nc = 0, Gtot = 0;
    const scat = p.enable_scattering && beta > 0;
    const phi = scat ? hgFovFraction(br.fov, p.dust_g_hg) : 0, w0 = p.dust_albedo;
    for (let j = 0; j < nl; j++) {
      if (h[j] === 0) continue;
      const tau = beta * d[j], tl = d[j] / C;
      _t[nc] = tl; _g[nc] = h[j] * Math.exp(-tau); Gtot += _g[nc]; nc++;
      if (scat) {
        const gs = h[j] * phi * (-Math.expm1(-w0 * tau)) * Math.exp(-(1 - w0) * tau);
        _t[nc] = tl + tl * (1 + tau) * (1 - p.dust_g_hg) * 0.5; _g[nc] = gs; Gtot += gs; nc++;
      }
    }
    // orden de las llegadas puntuales (pocas) y fusión con el histograma (ya ordenado en el tiempo)
    for (let i = 0; i < nc; i++) _ord[i] = i;
    for (let i = 1; i < nc; i++) { const v = _ord[i]; let k = i - 1; while (k >= 0 && _t[_ord[k]] > _t[v]) { _ord[k + 1] = _ord[k]; k--; } _ord[k + 1] = v; }
    let Hsum = 0;
    if (hist) for (let k = 0; k < nb; k++) { const g = hist[off + k]; if (g > 0) Hsum += g * Math.exp(-beta * C * (k + 0.5) * dt); }
    Gtot += Hsum;
    if (Gtot <= 0) return { As: 0, Al: 0, Ae: 0, v: noiseVar(p, 0, br.A).v, Pdc: 0 };
    // mediana ponderada de los tiempos de llegada
    const half = 0.5 * Gtot; let cum = 0, t0 = NaN, ip = 0, k = 0;
    while (isNaN(t0)) {
      const tp = ip < nc ? _t[_ord[ip]] : Infinity;
      let th = Infinity, gh = 0;
      if (hist) { while (k < nb && !(hist[off + k] > 0)) k++; if (k < nb) { th = k * dt; gh = hist[off + k] * Math.exp(-beta * C * (k + 0.5) * dt); } }
      if (tp === Infinity && th === Infinity) { t0 = 0; break; }
      if (tp <= th) { cum += _g[_ord[ip]]; if (cum >= half) t0 = tp; ip++; }
      else { cum += gh; if (cum >= half) t0 = th; k++; }
    }
    let S = 0, Il = 0, Ie = 0;
    const acc = (t, g) => { const dl = t - t0, w = Math.min(Math.abs(dl) / Tb, 1); S += g * (1 - w); if (dl > 0) Il += g * w; else if (dl < 0) Ie += g * w; };
    for (let i = 0; i < nc; i++) acc(_t[i], _g[i]);
    if (hist) for (let q = 0; q < nb; q++) { const g = hist[off + q]; if (g > 0) acc(q * dt, g * Math.exp(-beta * C * (q + 0.5) * dt)); }
    const eta = ledBwEff(p.Rb, D.f3db(p));
    Il += (1 - eta) * S; S *= eta;
    const amp = p.resp * p.mod_index * scale, Pdc = scale * Gtot;
    return { As: amp * S, Al: amp * Il, Ae: amp * Ie, v: noiseVar(p, Pdc, br.A).v, Pdc };
  }
  function eyeBer(As, Al, Ae, v) {
    const s = Math.sqrt(Math.max(v, 1e-300));
    return 0.25 * (qfunc((As + Al + Ae) / s) + qfunc((As + Al - Ae) / s) + qfunc((As - Al + Ae) / s) + qfunc((As - Al - Ae) / s));
  }
  /* Rendimiento en un punto: combina las ramas (MRC si ADR). */
  function pointPerf(p, brs, los, refl, i) {
    let as, al, ae, v;
    if (brs.length === 1) {
      const e = eyeTerms(p, brs[0], los[0].h, los[0].d, refl ? refl.hist[0] : null, i * (refl ? refl.nb : 0), refl ? refl.dt : 0, refl ? refl.nb : 0);
      as = e.As; al = e.Al; ae = e.Ae; v = Math.max(e.v, 1e-300);
    } else {
      let Ss = 0, L = 0, E = 0;
      for (let b = 0; b < brs.length; b++) {
        const e = eyeTerms(p, brs[b], los[b].h, los[b].d, refl ? refl.hist[b] : null, i * (refl ? refl.nb : 0), refl ? refl.dt : 0, refl ? refl.nb : 0);
        const vb = Math.max(e.v, 1e-300);
        Ss += e.As * e.As / vb; L += e.As * e.Al / vb; E += e.As * e.Ae / vb;
      }
      const nrm = Math.sqrt(Math.max(Ss, 1e-300));
      as = nrm; al = L / nrm; ae = E / nrm; v = 1;
    }
    const berPre = eyeBer(as, al, ae, v);
    const ber = p.use_fec ? rsPostFecBer(berPre) : berPre;
    const eye = Math.max(as - al - ae, 0);
    const snrDb = 10 * Math.log10(Math.max(eye * eye / v, 1e-30));
    const gp = goodputMbps(p.Rb, ber, D.codeRate(p), p.packet_bits);
    return { snrDb, berPre, ber, gp };
  }
  function evaluateVLC(p, grid, boxes, refl) {
    const brs = branches(p), n = grid.n;
    const snr = new Float32Array(n), ber = new Float64Array(n), gp = new Float32Array(n);
    for (let i = 0; i < n; i++) {
      const x = grid.pts[3 * i], y = grid.pts[3 * i + 1], z = grid.pts[3 * i + 2];
      const r = pointPerf(p, brs, losChannel(p, x, y, z, boxes, brs), refl, i);
      snr[i] = r.snrDb; ber[i] = r.ber; gp[i] = r.gp;
    }
    return { snr, ber, gp };
  }

  // ---------------------------------------------------------------- Wi-Fi 2,4 GHz (3GPP InF-SH + Middleton Class A)
  function wifiDefaults() {
    return { fc_GHz: 2.4, Pt_dBm: 20, G_ant_dBi: 2, BW: 20e6, Rb_target: 10e6, NF_dB: 7, ap_pos: [7.5, 7.5, 5.5],
      sigma_los: 4.3, sigma_nlos: 5.9, K_los_dB: 7, body_loss_dB: 4, emi_A: 0.1, emi_floor_dB: 14, emi_decay_exp: 2,
      use_fec: false, packet_bits: 1000, P_ap_W: 12, P_rx_wifi_W: 1.5, emi_scale_dB: 0 };
  }
  function emiAt(w, x, y, z) {
    let pr = Math.pow(10, w.emi_floor_dB / 10);
    for (const e of EMI_SOURCES) {
      const d = Math.max(Math.hypot(x - e.s[0], y - e.s[1], z - e.s[2]), 1);
      pr += Math.pow(10, e.s[3] / 10) * Math.pow(d, -w.emi_decay_exp);
    }
    return 10 * Math.log10(pr) + (w.emi_scale_dB || 0);
  }
  function middletonBer(g0, emiDb, A) {
    const inv = emiDb > -60 ? Math.pow(10, emiDb / 10) : 0;
    let ber = 0, pm = Math.exp(-A);
    for (let m = 0; m <= 20; m++) {
      if (m > 0) pm *= A / m;
      if (pm < 1e-15) break;
      ber += pm * qfunc(Math.sqrt(2 * g0 / (1 + (m / A) * inv)));
    }
    return ber;
  }
  /* Wi-Fi en un punto: disponibilidad sobre nFad muestras de desvanecimiento. */
  function wifiPoint(w, x, y, z, heavy, workers, sfNorm, r, nFad) {
    const ap = w.ap_pos, d3 = Math.max(Math.hypot(x - ap[0], y - ap[1], z - ap[2]), 1);
    const los = !(heavy.length && segBlocked(x, y, z, ap[0], ap[1], ap[2], heavy));
    const body = workers.length && segBlocked(x, y, z, ap[0], ap[1], ap[2], workers) ? w.body_loss_dB : 0;
    const plL = 31.84 + 21.5 * Math.log10(d3) + 19.0 * Math.log10(w.fc_GHz);
    const plN = Math.max(32.4 + 23.0 * Math.log10(d3) + 20.0 * Math.log10(w.fc_GHz), plL);
    const pl = los ? plL : plN, sf = (los ? w.sigma_los : w.sigma_nlos) * sfNorm;
    const prx = w.Pt_dBm + 2 * w.G_ant_dBi - pl - sf - body;
    const nTh = -174 + 10 * Math.log10(w.BW) + w.NF_dB, emi = emiAt(w, x, y, z);
    const K = Math.pow(10, w.K_los_dB / 10), s = Math.sqrt(K / (K + 1)), sg = Math.sqrt(1 / (2 * (K + 1)));
    let ok = 0, gpSum = 0;
    for (let f = 0; f < nFad; f++) {
      let pw;
      if (los) { const a = s + sg * r.normal(), b = sg * r.normal(); pw = a * a + b * b; }
      else { const a = r.normal() / Math.SQRT2, b = r.normal() / Math.SQRT2; pw = a * a + b * b; }
      const snrTh = prx + 10 * Math.log10(Math.max(pw, 1e-30)) - nTh;
      const g0 = Math.pow(10, snrTh / 10) * (w.BW / w.Rb_target);
      const bp = middletonBer(g0, emi, w.emi_A), bb = w.use_fec ? rsPostFecBer(bp) : bp;
      if (bb <= 1e-6) ok++;
      gpSum += goodputMbps(w.Rb_target, bb, w.use_fec ? 239 / 255 : 1, w.packet_bits);
    }
    return { avail: ok / nFad, gp: gpSum / nFad, los, emi, snrMeanDb: prx - (nTh + 10 * Math.log10(1 + Math.pow(10, emi / 10))) };
  }
  const wifiPcomm = (w) => w.P_ap_W + w.P_rx_wifi_W;

  // ---------------------------------------------------------------- iluminancia (EN 12464-1)
  function illuminance(p, grid, withReflections = true) {
    const m = D.m(p), I0 = p.P_opt * p.LER * (m + 1) / (2 * Math.PI), E = new Float64Array(grid.n);
    for (let i = 0; i < grid.n; i++) {
      const x = grid.pts[3 * i], y = grid.pts[3 * i + 1], z = grid.pts[3 * i + 2];
      for (const [lx, ly] of p.led_xy) { const vz = p.z_led - z, d = Math.hypot(lx - x, ly - y, vz), c = vz / d; if (c > 0) E[i] += I0 * Math.pow(c, m) * c / (d * d); }
    }
    if (withReflections) {
      const W = walls(p), Ew = new Float64Array(W.c.length);
      W.c.forEach((c, w) => { for (const [lx, ly] of p.led_xy) { const vx = c[0] - lx, vy = c[1] - ly, vz = c[2] - p.z_led, d = Math.hypot(vx, vy, vz);
        const cp = -vz / d, ca = -(vx * W.n[w][0] + vy * W.n[w][1] + vz * W.n[w][2]) / d; if (cp > 0 && ca > 0) Ew[w] += I0 * Math.pow(cp, m) * ca / (d * d); } });
      for (let i = 0; i < grid.n; i++) {
        const x = grid.pts[3 * i], y = grid.pts[3 * i + 1], z = grid.pts[3 * i + 2];
        W.c.forEach((c, w) => { const vx = x - c[0], vy = y - c[1], vz = z - c[2], d = Math.hypot(vx, vy, vz);
          const cb = (vx * W.n[w][0] + vy * W.n[w][1] + vz * W.n[w][2]) / d, cpsi = -vz / d;
          if (cb > 0 && cpsi > 0) E[i] += p.rho_wall * Ew[w] * W.dA / Math.PI * cb * cpsi / (d * d); });
      }
    }
    let s = 0, mn = Infinity; for (const e of E) { s += e; mn = Math.min(mn, e); }
    return { E, mean: s / grid.n, u0: mn / (s / grid.n) };
  }

  return { qfunc, rsPostFecBer, goodputMbps, hgFovFraction, ledBwEff, defaults, optimizedGrid, D, STATIC, STATIC_BOXES,
    EMI_SOURCES, rng, randomDynamic, segBlocked, insideAny, receiverGrid, walls, branches, losChannel, reflections,
    eyeTerms, eyeBer, pointPerf, evaluateVLC, wifiDefaults, wifiPoint, wifiPcomm, emiAt, middletonBer, illuminance,
    GRID3, GRID4 };
})();
if (typeof module !== "undefined") module.exports = VLC;
