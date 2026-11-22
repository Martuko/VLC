# -*- coding: utf-8 -*-
"""
vlc_core.py
===========
Núcleo del simulador de un sistema VLC (Visible Light Communication) en una
nave industrial hostil.

Modelos implementados
---------------------
1. Canal óptico LOS Lambertiano + reflexiones de primer orden en las paredes
   (Komine & Nakagawa, IEEE Trans. Consumer Electronics, 2004;
    Ghassemlooy et al., "Optical Wireless Communications", CRC Press).
2. Atenuación por polvo en suspensión (ley de Beer-Lambert) con coeficiente de
   extinción obtenido de la visibilidad (relación de Koschmieder: beta = 3.912 / V).
3. Suciedad depositada sobre las ópticas (cubierta de la luminaria y lente del
   receptor): transmitancia que decae exponencialmente con los días sin limpieza.
4. Temperatura: caída del flujo del LED con la temperatura, aumento del ruido
   térmico del preamplificador y de la corriente oscura del fotodiodo.
5. Obstáculos fijos (máquinas, estanterías) y móviles (operarios, grúas horquilla)
   que bloquean trayectorias (prueba segmento-caja por el método de "slabs").
6. Ruido shot + térmico (modelo de Komine), interferencia intersimbólica (ISI)
   por la dispersión temporal del canal, BER de OOK-NRZ.
7. Referencia RF: Wi-Fi 2,4 GHz con el modelo de pérdidas 3GPP TR 38.901
   "Indoor Factory" (InF-SH), sombreamiento log-normal, desvanecimiento
   Rice/Rayleigh e interferencia electromagnética (EMI) industrial.

Todas las unidades están en el Sistema Internacional (m, W, A, s, K) salvo que
se indique lo contrario.
"""
from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field, replace
from scipy.special import erfc

# ---------------------------------------------------------------------------
# Constantes físicas
# ---------------------------------------------------------------------------
Q_E = 1.602176634e-19      # carga del electrón [C]
K_B = 1.380649e-23         # constante de Boltzmann [J/K]
C_LUZ = 2.99792458e8       # velocidad de la luz [m/s]


def qfunc(x):
    """Función Q gaussiana."""
    return 0.5 * erfc(np.asarray(x) / np.sqrt(2.0))


# ---------------------------------------------------------------------------
# Parámetros del sistema
# ---------------------------------------------------------------------------
@dataclass
class Params:
    # --- Recinto (nave industrial) ---
    L: float = 15.0              # largo  (x) [m]
    W: float = 15.0              # ancho  (y) [m]
    H: float = 6.0               # alto   (z) [m]
    rho_wall: float = 0.4        # reflectividad de las paredes (hormigón/chapa sucia)
    dA_wall: float = 0.5         # tamaño del elemento de pared para reflexiones [m]

    # --- Transmisores (luminarias LED tipo "high-bay") ---
    led_xy: tuple = ((2.5, 2.5), (2.5, 7.5), (2.5, 12.5),
                     (7.5, 2.5), (7.5, 7.5), (7.5, 12.5),
                     (12.5, 2.5), (12.5, 7.5), (12.5, 12.5))
    z_led: float = 5.8           # altura de montaje [m]
    P_opt: float = 45.0          # potencia óptica media por luminaria [W]
    LER: float = 300.0           # eficacia luminosa de la radiación del LED blanco [lm/W]
    phi_half_deg: float = 60.0   # semiángulo de media potencia [°]
    mod_index: float = 0.2       # índice de modulación (amplitud AC / potencia media)

    # --- Receptor (fotodiodo PIN + concentrador) ---
    z_rx: float = 1.0            # altura del plano de recepción (AGV / sensor) [m]
    A_pd: float = 1.0e-4         # área del fotodiodo [m^2] (1 cm^2)
    fov_deg: float = 60.0        # campo de visión (FOV) [°]
    n_conc: float = 1.5          # índice de refracción del concentrador
    Ts_filter: float = 1.0       # transmitancia del filtro óptico
    resp: float = 0.53           # responsividad [A/W]

    # --- Ruido (parámetros de Komine & Nakagawa 2004) ---
    I_bg: float = 1000e-6        # corriente por luz ambiental de fondo [A]
    I2: float = 0.562            # factor de ancho de banda de ruido
    I3: float = 0.0868
    G_ol: float = 10.0           # ganancia de lazo abierto del preamplificador
    eta_cap: float = 112e-12 / 1e-4   # capacitancia por unidad de área [F/m^2] (112 pF/cm^2)
    Gamma: float = 1.5           # factor de ruido del canal del FET
    g_m: float = 30e-3           # transconductancia del FET [S]
    I_dark0: float = 10e-9       # corriente oscura a 25 °C [A]

    # --- Enlace ---
    Rb: float = 50e6             # tasa de bits OOK-NRZ [bit/s]  (B = Rb)

    # --- Condiciones hostiles ---
    T_amb: float = 25.0          # temperatura ambiente [°C]
    k_T_led: float = 0.004       # caída del flujo del LED [1/°C] (≈ -0,4 %/°C)
    visibility: float = np.inf   # visibilidad óptica por polvo [m] (inf = aire limpio)
    days_dirty: float = 0.0      # días sin limpiar las ópticas
    tau_rx_days: float = 30.0    # constante de ensuciamiento del lente del receptor (mira hacia arriba)
    tau_tx_days: float = 120.0   # constante de ensuciamiento de la cubierta de la luminaria (mira hacia abajo)

    # ----- Propiedades derivadas -----
    @property
    def m_lamb(self):
        """Orden Lambertiano m = -ln2 / ln(cos(Phi_1/2))."""
        return -np.log(2.0) / np.log(np.cos(np.radians(self.phi_half_deg)))

    @property
    def fov(self):
        return np.radians(self.fov_deg)

    @property
    def g_conc(self):
        """Ganancia del concentrador no formador de imagen."""
        return self.n_conc ** 2 / np.sin(self.fov) ** 2

    @property
    def beta_ext(self):
        """Coeficiente de extinción por polvo [1/m] (Koschmieder, 550 nm)."""
        return 0.0 if np.isinf(self.visibility) else 3.912 / self.visibility

    @property
    def T_soil(self):
        """Transmitancia total por suciedad depositada (luminaria × receptor)."""
        return np.exp(-self.days_dirty / self.tau_rx_days) * \
               np.exp(-self.days_dirty / self.tau_tx_days)

    # Protección térmica del driver (opcional): sobre T_derate_C la corriente se
    # reduce linealmente con pendiente derate_slope [1/°C]. None = sin protección.
    T_derate_C: float = None
    derate_slope: float = 0.02

    @property
    def led_derating(self):
        """Factor de potencia óptica relativa por temperatura."""
        f = 1.0 - self.k_T_led * (self.T_amb - 25.0)
        if self.T_derate_C is not None and self.T_amb > self.T_derate_C:
            f *= 1.0 - self.derate_slope * (self.T_amb - self.T_derate_C)
        return max(0.0, f)

    @property
    def leds(self):
        xy = np.array(self.led_xy, dtype=float)
        return np.column_stack([xy, np.full(len(xy), self.z_led)])


# ---------------------------------------------------------------------------
# Obstáculos: cajas alineadas a los ejes [xmin, xmax, ymin, ymax, zmax] (zmin = 0)
# ---------------------------------------------------------------------------
STATIC_OBSTACLES = {
    # nombre: (xmin, xmax, ymin, ymax, altura)
    "Estantería A":        (1.0, 7.0, 13.2, 14.2, 5.0),
    "Estantería B":        (8.0, 14.0, 13.2, 14.2, 5.0),
    "Prensa":              (1.5, 3.5, 1.5, 3.5, 4.0),
    "CNC 1":               (5.5, 8.5, 2.5, 4.5, 2.5),
    "CNC 2":               (10.5, 13.5, 1.5, 4.0, 2.5),
    "Horno":               (10.0, 13.0, 7.0, 10.0, 3.5),
    "Cinta transportadora": (1.0, 8.5, 8.0, 8.8, 1.3),
}


def static_boxes():
    return np.array(list(STATIC_OBSTACLES.values()), dtype=float)


def random_dynamic_boxes(rng, n_workers, n_forklifts, p: Params, static=None):
    """Genera operarios (0,5×0,5×1,8 m) y grúas horquilla (1,2×2,4×2,3 m)
    en posiciones aleatorias fuera de las máquinas."""
    boxes = []
    static = static_boxes() if static is None else static

    def free(x0, x1, y0, y1):
        for b in static:
            if not (x1 < b[0] or x0 > b[1] or y1 < b[2] or y0 > b[3]):
                return False
        return True

    for (n, sx, sy, h) in ((n_workers, 0.5, 0.5, 1.8), (n_forklifts, 1.2, 2.4, 2.3)):
        k = 0
        while k < n:
            # orientación aleatoria de la grúa horquilla
            if rng.random() < 0.5:
                sx_, sy_ = sx, sy
            else:
                sx_, sy_ = sy, sx
            cx = rng.uniform(sx_ / 2, p.L - sx_ / 2)
            cy = rng.uniform(sy_ / 2, p.W - sy_ / 2)
            b = (cx - sx_ / 2, cx + sx_ / 2, cy - sy_ / 2, cy + sy_ / 2, h)
            if free(*b[:4]):
                boxes.append(b)
                k += 1
    return np.array(boxes, dtype=float).reshape(-1, 5)


def segment_blocked(P0, P1, boxes, ignore_if_inside=True):
    """
    Determina si los segmentos P0->P1 intersectan alguna caja.
    P0, P1 : (N,3) ; boxes : (K,5)
    Método de "slabs" (Kay-Kajiya). Si el extremo P0 (receptor) está dentro de
    la huella de la caja, esa caja se ignora (el operario/la máquina porta el equipo).
    Devuelve un arreglo booleano (N,).
    """
    P0 = np.atleast_2d(P0)
    P1 = np.atleast_2d(P1)
    N = P0.shape[0]
    blocked = np.zeros(N, dtype=bool)
    if boxes is None or len(boxes) == 0:
        return blocked
    D = P1 - P0
    eps = 1e-12
    invD = 1.0 / np.where(np.abs(D) < eps, np.where(D >= 0, eps, -eps), D)
    for b in boxes:
        lo = np.array([b[0], b[2], 0.0])
        hi = np.array([b[1], b[3], b[4]])
        t1 = (lo - P0) * invD
        t2 = (hi - P0) * invD
        tmin = np.max(np.minimum(t1, t2), axis=1)
        tmax = np.min(np.maximum(t1, t2), axis=1)
        hit = (tmax >= np.maximum(tmin, 0.0)) & (tmin <= 1.0) & (tmax > 1e-6)
        if ignore_if_inside:
            inside = (P0[:, 0] >= b[0]) & (P0[:, 0] <= b[1]) & \
                     (P0[:, 1] >= b[2]) & (P0[:, 1] <= b[3])
            hit &= ~inside
        blocked |= hit
    return blocked


def inside_any(xy, boxes):
    """Máscara de puntos (x,y) que caen dentro de la huella de alguna caja."""
    m = np.zeros(len(xy), dtype=bool)
    for b in boxes:
        m |= (xy[:, 0] >= b[0]) & (xy[:, 0] <= b[1]) & (xy[:, 1] >= b[2]) & (xy[:, 1] <= b[3])
    return m


# ---------------------------------------------------------------------------
# Geometría auxiliar
# ---------------------------------------------------------------------------
def receiver_grid(p: Params, step=0.25):
    xs = np.arange(step / 2, p.L, step)
    ys = np.arange(step / 2, p.W, step)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, p.z_rx)])
    return pts, xs, ys


def wall_elements(p: Params):
    """Centros, normales (hacia el interior) y área de los elementos de pared."""
    d = p.dA_wall
    cs, ns = [], []
    zs = np.arange(d / 2, p.H, d)
    xs = np.arange(d / 2, p.L, d)
    ys = np.arange(d / 2, p.W, d)
    for x in xs:
        for z in zs:
            cs.append((x, 0.0, z)); ns.append((0, 1, 0))      # pared y = 0
            cs.append((x, p.W, z)); ns.append((0, -1, 0))     # pared y = W
    for y in ys:
        for z in zs:
            cs.append((0.0, y, z)); ns.append((1, 0, 0))      # pared x = 0
            cs.append((p.L, y, z)); ns.append((-1, 0, 0))     # pared x = L
    return np.array(cs), np.array(ns, dtype=float), d * d


# ---------------------------------------------------------------------------
# Canal óptico
# ---------------------------------------------------------------------------
def illuminance(p: Params, rx):
    """Iluminancia horizontal directa [lx] en los puntos rx."""
    m = p.m_lamb
    I0 = p.P_opt * p.LER * (m + 1) / (2 * np.pi)       # intensidad axial [cd]
    E = np.zeros(len(rx))
    for led in p.leds:
        v = led - rx
        d = np.linalg.norm(v, axis=1)
        cos = v[:, 2] / d
        E += I0 * cos ** m * cos / d ** 2
    return E


@dataclass
class ChannelResult:
    """Resultado del canal (ganancias sin incluir polvo/suciedad/temperatura)."""
    h_los: np.ndarray        # (Nrx, Nled) ganancia DC LOS (0 si bloqueado o fuera de FOV)
    d_los: np.ndarray        # (Nrx, Nled) distancia LOS [m]
    los_ok: np.ndarray       # (Nrx, Nled) LOS geométricamente posible (antes de bloqueo dinámico)
    hist: dict = field(default_factory=dict)   # histograma de reflexiones por beta_ext
    t_bins: np.ndarray = None


def los_channel(p: Params, rx, boxes=None):
    """Ganancia DC de línea de vista para cada par receptor-luminaria."""
    leds = p.leds
    m, A, g = p.m_lamb, p.A_pd, p.g_conc
    Nrx, Nl = len(rx), len(leds)
    h = np.zeros((Nrx, Nl))
    dist = np.zeros((Nrx, Nl))
    for j, led in enumerate(leds):
        v = led - rx
        d = np.linalg.norm(v, axis=1)
        cos_phi = v[:, 2] / d          # LED apunta hacia abajo
        cos_psi = v[:, 2] / d          # receptor apunta hacia arriba
        psi = np.arccos(np.clip(cos_psi, -1, 1))
        hj = (m + 1) * A / (2 * np.pi * d ** 2) * cos_phi ** m * p.Ts_filter * g * cos_psi
        hj[psi > p.fov] = 0.0
        if boxes is not None and len(boxes):
            blk = segment_blocked(rx, np.repeat(led[None, :], Nrx, 0), boxes)
            hj[blk] = 0.0
        h[:, j] = hj
        dist[:, j] = d
    return h, dist


def reflection_histogram(p: Params, rx, boxes=None, betas=(0.0,), dt=0.5e-9, t_max=300e-9,
                         chunk=150, verbose=False):
    """
    Respuesta impulsiva de las reflexiones de primer orden, discretizada en
    intervalos de dt. Se calcula para varios coeficientes de extinción beta
    (polvo), porque la atenuación depende de la longitud de cada trayecto.
    Devuelve dict beta -> (Nrx, Nbins) y el vector de tiempos.
    """
    leds = p.leds
    m, A, g, rho = p.m_lamb, p.A_pd, p.g_conc, p.rho_wall
    wc, wn, dA = wall_elements(p)
    nb = int(np.ceil(t_max / dt))
    t_bins = np.arange(nb) * dt
    out = {b: np.zeros((len(rx), nb)) for b in betas}

    # --- tramo luminaria -> pared (no depende del receptor) ---
    seg1 = []
    for led in leds:
        v1 = wc - led                       # desde el LED hacia la pared
        d1 = np.linalg.norm(v1, axis=1)
        cos_phi = -v1[:, 2] / d1            # ángulo respecto del eje del LED (0,0,-1)
        cos_alpha = np.einsum("ij,ij->i", -v1, wn) / d1   # incidencia en la pared
        ok = (cos_phi > 0) & (cos_alpha > 0)
        if boxes is not None and len(boxes):
            ok &= ~segment_blocked(wc, np.repeat(led[None, :], len(wc), 0), boxes, ignore_if_inside=False)
        a1 = np.where(ok, (m + 1) / (2 * np.pi) * cos_phi ** m * cos_alpha / d1 ** 2, 0.0)
        seg1.append((a1, d1))

    # --- tramo pared -> receptor ---
    for s in range(0, len(rx), chunk):
        r = rx[s:s + chunk]                                  # (n,3)
        v2 = r[:, None, :] - wc[None, :, :]                  # pared -> rx  (n,Nw,3)
        d2 = np.linalg.norm(v2, axis=2)
        cos_beta = np.einsum("nwk,wk->nw", v2, wn) / d2      # salida desde la pared
        cos_psi = -v2[:, :, 2] / d2                          # llegada al rx (mira hacia arriba)
        psi = np.arccos(np.clip(cos_psi, -1, 1))
        ok = (cos_beta > 0) & (cos_psi > 0) & (psi <= p.fov)
        if boxes is not None and len(boxes):
            P0 = np.repeat(r, len(wc), 0)
            P1 = np.tile(wc, (len(r), 1))
            blk = segment_blocked(P0, P1, boxes).reshape(len(r), len(wc))
            ok &= ~blk
        a2 = np.where(ok, rho * dA * cos_beta / (np.pi * d2 ** 2) * A * p.Ts_filter * g * cos_psi, 0.0)
        for (a1, d1) in seg1:
            gain = a1[None, :] * a2                          # (n,Nw)
            dtot = d1[None, :] + d2
            idx = np.minimum((dtot / C_LUZ / dt).astype(int), nb - 1)
            for b in betas:
                gb = gain * np.exp(-b * dtot)
                for i in range(len(r)):
                    out[b][s + i] += np.bincount(idx[i], weights=gb[i], minlength=nb)
        if verbose:
            print(f"   reflexiones: {min(s + chunk, len(rx))}/{len(rx)} receptores", flush=True)
    return out, t_bins


# ---------------------------------------------------------------------------
# Rendimiento: potencia, SNR, BER
# ---------------------------------------------------------------------------
def noise_variance(p: Params, P_dc_rx):
    """
    Varianza total de ruido [A^2] (Komine 2004) con dependencia de temperatura.
    P_dc_rx : potencia óptica media total recibida [W] (genera ruido shot).
    """
    B = p.Rb
    T_K = p.T_amb + 273.15
    I_dark = p.I_dark0 * 2.0 ** ((p.T_amb - 25.0) / 10.0)
    s_shot = 2 * Q_E * p.resp * P_dc_rx * B + 2 * Q_E * p.I_bg * p.I2 * B + 2 * Q_E * I_dark * B
    C = p.eta_cap * p.A_pd
    s_th = (8 * np.pi * K_B * T_K / p.G_ol) * C * p.I2 * B ** 2 + \
           (16 * np.pi ** 2 * K_B * T_K * p.Gamma / p.g_m) * C ** 2 * p.I3 * B ** 3
    return s_shot, s_th


def performance(p: Params, h_los, d_los, refl_hist=None, t_bins=None):
    """
    Calcula potencia recibida, SNR y BER por receptor.
    h_los, d_los : (Nrx, Nled)   refl_hist : (Nrx, Nbins) ya atenuado por polvo.
    """
    beta = p.beta_ext
    scale = p.P_opt * p.led_derating * p.T_soil        # potencia media efectiva por luminaria
    att_los = np.exp(-beta * d_los)
    g_los = h_los * att_los                              # ganancia LOS con polvo
    t_los = d_los / C_LUZ

    Nrx = h_los.shape[0]
    Tb = 1.0 / p.Rb
    G_sig = np.zeros(Nrx)
    G_isi = np.zeros(Nrx)

    # instante de primera llegada
    has_los = g_los.max(axis=1) > 0
    t0 = np.where(has_los, np.where(g_los > 0, t_los, np.inf).min(axis=1), np.inf)
    if refl_hist is not None:
        first_refl = np.argmax(refl_hist > 0, axis=1)
        any_refl = refl_hist.max(axis=1) > 0
        t0 = np.where(np.isinf(t0) & any_refl, t_bins[first_refl], t0)

    # componentes LOS: dentro o fuera del periodo de bit
    in_bit = (t_los - t0[:, None]) <= Tb
    G_sig += np.sum(np.where(in_bit, g_los, 0.0), axis=1)
    G_isi += np.sum(np.where(in_bit, 0.0, g_los), axis=1)

    if refl_hist is not None:
        rel = t_bins[None, :] - t0[:, None]
        mask = (rel >= -1e-12) & (rel <= Tb)
        G_sig += np.sum(np.where(mask, refl_hist, 0.0), axis=1)
        G_isi += np.sum(np.where(mask, 0.0, refl_hist), axis=1)

    G_tot = G_sig + G_isi
    P_dc = scale * G_tot                                 # potencia óptica media recibida
    I_sig = p.resp * p.mod_index * scale * G_sig          # amplitud de la señal [A]
    I_isi = p.resp * p.mod_index * scale * G_isi
    s_shot, s_th = noise_variance(p, P_dc)
    snr = I_sig ** 2 / (s_shot + s_th + I_isi ** 2)
    ber = qfunc(np.sqrt(snr))
    return dict(P_dc=P_dc, snr=snr, snr_db=10 * np.log10(np.maximum(snr, 1e-30)),
                ber=ber, s_shot=s_shot, s_th=s_th, G_los=g_los.sum(1), G_tot=G_tot)


def rms_delay_spread(h_los, d_los, refl_hist, t_bins):
    """Dispersión de retardo RMS [s] por receptor."""
    t_los = d_los / C_LUZ
    num1 = (h_los * t_los).sum(1) + (refl_hist * t_bins[None, :]).sum(1)
    num2 = (h_los * t_los ** 2).sum(1) + (refl_hist * t_bins[None, :] ** 2).sum(1)
    den = h_los.sum(1) + refl_hist.sum(1)
    mu = num1 / den
    return np.sqrt(np.maximum(num2 / den - mu ** 2, 0))


# ---------------------------------------------------------------------------
# Referencia RF: Wi-Fi 2,4 GHz en fábrica (3GPP TR 38.901 InF-SH)
# ---------------------------------------------------------------------------
@dataclass
class WifiParams:
    fc_GHz: float = 2.4
    Pt_dBm: float = 20.0
    G_ant_dBi: float = 2.0          # ganancia de cada antena
    BW: float = 20e6
    NF_dB: float = 7.0
    ap_pos: tuple = (7.5, 7.5, 5.5)
    sigma_los: float = 4.3
    sigma_nlos: float = 5.9
    K_los_dB: float = 7.0            # factor de Rice en LOS (supuesto)
    emi_dB: float = 0.0              # aumento del piso de ruido por EMI [dB sobre ruido térmico]


def wifi_pathloss(w: WifiParams, d3d, los):
    pl_los = 31.84 + 21.50 * np.log10(d3d) + 19.00 * np.log10(w.fc_GHz)
    pl_sh = 32.4 + 23.0 * np.log10(d3d) + 20.0 * np.log10(w.fc_GHz)
    pl_nlos = np.maximum(pl_sh, pl_los)
    return np.where(los, pl_los, pl_nlos)


def wifi_noise_dBm(w: WifiParams):
    n_th = -174 + 10 * np.log10(w.BW) + w.NF_dB
    return n_th + 10 * np.log10(1 + 10 ** (w.emi_dB / 10))


def wifi_fading_power(rng, los, n):
    """Ganancia de potencia |h|^2 con media 1: Rice (LOS, K) o Rayleigh (NLOS)."""
    K = 10 ** (7.0 / 10)
    out = np.empty((len(los), n))
    for i, l in enumerate(los):
        if l:
            s = np.sqrt(K / (K + 1))
            sig = np.sqrt(1 / (2 * (K + 1)))
            hr = s + sig * rng.standard_normal(n)
            hi = sig * rng.standard_normal(n)
        else:
            hr = rng.standard_normal(n) / np.sqrt(2)
            hi = rng.standard_normal(n) / np.sqrt(2)
        out[i] = hr ** 2 + hi ** 2
    return out


# Umbrales de SNR para BER = 1e-6
SNR_TH_OOK_DB = 10 * np.log10(4.7534 ** 2)          # OOK: BER = Q(sqrt(SNR))   -> 13,5 dB
SNR_TH_BPSK_DB = 10 * np.log10(4.7534 ** 2 / 2)     # BPSK: BER = Q(sqrt(2 SNR)) -> 10,5 dB
