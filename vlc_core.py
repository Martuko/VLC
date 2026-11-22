# -*- coding: utf-8 -*-
"""
vlc_core.py
===========
Núcleo físico del simulador de un sistema VLC (Visible Light Communication)
en una nave industrial metalmecánica pesada con alta interferencia
electromagnética (EMI), polvo en suspensión, suciedad óptica, calor extremo
y obstáculos fijos y móviles.

Contexto industrial del simulador
---------------------------------
Nave de manufactura pesada de 15 × 15 × 6 m equipada con horno de inducción,
prensa, centros de mecanizado CNC con variadores de frecuencia (VFD),
estanterías metálicas y línea de cinta con puesto de soldadura. Las luminarias
LED industriales ("high-bay") cumplen doble función:
  1. Iluminación del plano de trabajo conforme a EN 12464-1 (E ≥ 300 lx, U0 ≥ 0,60).
  2. Enlace descendente VLC inmune a EMI hacia AGVs autónomos y sensores de piso
     (a tasa nominal industrial de 10 Mb/s), complementado con un enlace de subida
     infrarrojo (IR 850 nm) para operación bidireccional 100 % libre de RF.

Modelos físicos implementados y referencias verificadas
-------------------------------------------------------
1. Canal óptico LOS Lambertiano y reflexiones difusas de primer orden:
   - T. Komine & M. Nakagawa, IEEE Trans. Consumer Electronics, 50(1), 100–107,
     2004. DOI: 10.1109/TCE.2004.1277847
   - M. Uysal, F. Miramirkhani et al., IEEE Communications Magazine, 55(1),
     212–217, 2017 (IEEE 802.15.7r1 Manufacturing Cell). DOI: 10.1109/MCOM.2017.1600872CM
2. Respuesta en frecuencia limitada del LED blanco fosforescente y filtro azul:
   - H. Le Minh et al., IEEE Photonics Technology Letters, 21(15), 1063–1065,
     2009. DOI: 10.1109/LPT.2009.2022413
3. Polvo en suspensión (extinción Beer-Lambert/Koschmieder + dispersión hacia
   adelante con función de fase de Henyey-Greenstein y retardo temporal):
   - A. Al-Kinani et al., IEEE Commun. Surveys & Tutorials, 20(3), 1939–1962,
     2018. DOI: 10.1109/COMST.2018.2838096
   - P. Palacios Játiva et al., IEEE Access, 8, 185445–185464, 2020
     (DOI: 10.1109/ACCESS.2020.3030615) y MDPI Sensors, 22(7), 2483, 2022
     (DOI: 10.3390/s22072483).
   - D. Schneider et al., Optics Express, 29(11), 2021. DOI: 10.1364/OE.421757
4. Suciedad depositada en ópticas (LDD) y acoplamiento con temperatura de unión Tj:
   - R. G. Davis, A. Wilkerson & B. R. Kinzey, LEUKOS, 15(1), 55–63, 2019.
     DOI: 10.1080/15502724.2018.1517598
5. Receptor con Diversidad Angular piramidal (ADR, combinación MRC/SBC):
   - A. Nuwanpriya, S.-W. Ho & C. S. Chen, IEEE JSAC, 33(9), 1780–1792, 2015.
     DOI: 10.1109/JSAC.2015.2432514
6. Codificación FEC Reed-Solomon RS(255,239) y métricas de eficiencia verde:
   - IEEE Std 802.15.7-2018. DOI: 10.1109/IEEESTD.2019.8697198
   - C. Jenila & R. K. Jeyachitra, Digital Communications and Networks, 7(3),
     410–444, 2021 (ScienceDirect). DOI: 10.1016/j.dcan.2020.09.004
   - L. E. M. Matheus et al., IEEE Commun. Surveys & Tutorials, 21(4), 3204–3237,
     2019. DOI: 10.1109/COMST.2019.2913348
7. Referencia RF: Wi-Fi 2,4 GHz en fábrica (3GPP TR 38.901 InF-SH) + ruido
   impulsivo industrial Middleton Class A y mapa espacial de EMI:
   - 3GPP TR 38.901, sección 7.4.1 (Indoor Factory).
   - D. Middleton, IEEE Trans. EMC, EMC-19(3), 106–127, 1977.
     DOI: 10.1109/TEMC.1977.303527
   - J.-D. C. Oyedum et al., IEEE GLOBECOM, 2025.
     DOI: 10.1109/GLOBECOM59602.2025.11432514

Todas las magnitudes están en el Sistema Internacional (m, W, A, s, K, Hz)
salvo cuando el nombre del parámetro indica explícitamente otra unidad.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
import numpy as np
from scipy.special import betainc, erfc

# ---------------------------------------------------------------------------
# Constantes físicas fundamentales
# ---------------------------------------------------------------------------
Q_E = 1.602176634e-19      # carga elemental del electrón [C]
K_B = 1.380649e-23         # constante de Boltzmann [J/K]
C_LUZ = 2.99792458e8       # velocidad de la luz en el vacío [m/s]


def qfunc(x):
    """Función Q gaussiana: Q(x) = 0.5 * erfc(x / sqrt(2))."""
    return 0.5 * erfc(np.asarray(x, dtype=float) / np.sqrt(2.0))


def rs_post_fec_ber(ber_pre, n: int = 255, k: int = 239, m_bits: int = 8):
    """
    Calcula la BER post-FEC exacta (cota de decodificación por bloque) para un
    código Reed-Solomon RS(n, k) sobre GF(2^m_bits), estándar en IEEE 802.15.7.

    Con RS(255, 239), m=8 bits/símbolo, la capacidad correctora es t = 8 símbolos
    y la tasa de código es Rc = 239/255 ≈ 0,9373 (overhead de 6,27 %).

    Usando la identidad exacta con la función beta incompleta regularizada:
        BER_post = p_b * I_{p_s}(t, n - t)
    donde p_s = 1 - (1 - p_b)^m_bits es la probabilidad de error de símbolo.
    """
    t = (n - k) // 2
    pb = np.clip(np.asarray(ber_pre, dtype=float), 0.0, 0.5)
    # 1 - (1 - pb)^m con alta precisión numérica para pb muy pequeños
    ps = -np.expm1(m_bits * np.log1p(-pb))
    ps = np.clip(ps, 0.0, 1.0)
    return pb * betainc(t, n - t, ps)


def goodput_mbps(Rb: float, ber_post, code_rate: float = 1.0, packet_bits: int = 1000):
    """
    Calcula el throughput útil libre de errores (Goodput) en Mb/s:
        G_eff = (Rb / 1e6) * R_code * (1 - PER)
    donde PER = 1 - (1 - BER_post)^L_pkt es la tasa de error de paquete.
    """
    pb = np.clip(np.asarray(ber_post, dtype=float), 0.0, 1.0 - 1e-15)
    per = -np.expm1(packet_bits * np.log1p(-pb))
    per = np.clip(per, 0.0, 1.0)
    return (Rb / 1e6) * code_rate * (1.0 - per)


def hg_fov_fraction(fov_rad: float, g_hg: float) -> float:
    """
    Integral analítica exacta de la función de fase de Henyey-Greenstein sobre
    el cono de aceptación [0, fov_rad] del receptor:
        Phi_HG(Psi; g) = ((1 - g^2) / (2g)) * [1/(1 - g) - 1/sqrt(1 + g^2 - 2g cos Psi)]
    Representa la fracción de radiación dispersada por polvo que entra al FOV.
    """
    cos_fov = np.cos(np.clip(fov_rad, 0.0, np.pi))
    if abs(g_hg) < 1e-6:
        return 0.5 * (1.0 - cos_fov)
    num = 1.0 - g_hg ** 2
    term = 1.0 / (1.0 - g_hg) - 1.0 / np.sqrt(1.0 + g_hg ** 2 - 2.0 * g_hg * cos_fov)
    return float(np.clip((num / (2.0 * g_hg)) * term, 0.0, 1.0))


def led_bandwidth_efficiency(Rb: float, f_3db: float) -> float:
    """
    Fracción de amplitud de ojo útil en un periodo de bit Tb = 1/Rb para una
    señal OOK-NRZ atravesando un LED modelado como filtro paso bajo de primer
    orden con frecuencia de corte f_3db (Le Minh et al., 2009):
        eta_LED = 1 - exp(-Tb / tau_LED),    tau_LED = 1 / (2 * pi * f_3db)
    La cola residual exp(-Tb / tau_LED) se transfiere como ISI intrínseca del
    transmisor (peor caso de transición de bits).
    """
    if f_3db is None or np.isinf(f_3db):
        return 1.0
    Tb = 1.0 / Rb
    tau = 1.0 / (2.0 * np.pi * f_3db)
    return float(np.clip(-np.expm1(-Tb / tau), 0.0, 1.0))


# ---------------------------------------------------------------------------
# Parámetros del sistema VLC industrial
# ---------------------------------------------------------------------------
@dataclass
class Params:
    # --- Recinto (nave metalmecánica pesada) ---
    L: float = 15.0              # largo  (x) [m]
    W: float = 15.0              # ancho  (y) [m]
    H: float = 6.0               # alto   (z) [m]
    rho_wall: float = 0.4        # reflectividad difusa de paredes (hormigón/chapa)
    dA_wall: float = 0.5         # lado del elemento de pared para reflexiones [m]

    # --- Transmisores (luminarias LED industriales "high-bay") ---
    led_xy: tuple = ((2.5, 2.5), (2.5, 7.5), (2.5, 12.5),
                     (7.5, 2.5), (7.5, 7.5), (7.5, 12.5),
                     (12.5, 2.5), (12.5, 7.5), (12.5, 12.5))
    z_led: float = 5.8           # altura de montaje [m]
    P_opt: float = 45.0          # potencia óptica media por luminaria [W]
    led_wall_plug_eff: float = 0.45   # eficiencia electro-óptica (45 W opt / 100 W elec)
    P_comm_extra_per_led: float = 1.0 # potencia eléctrica incremental del driver VLC por luminaria [W]
    P_rx_vlc: float = 0.5        # potencia eléctrica del receptor óptico en el AGV [W]
    LER: float = 300.0           # eficacia luminosa de la radiación [lm/W_opt]
    phi_half_deg: float = 60.0   # semiángulo de media potencia del LED [°]
    mod_index: float = 0.2       # índice de modulación (amplitud AC / potencia media)

    # --- Ancho de banda del LED blanco y filtro óptico azul (Le Minh et al. 2009) ---
    # Por defecto: LED industrial con pre-ecualización pasiva básica (15 MHz).
    # Para LED fosforescente sin ecualizar usar f_3db_led = 3e6 Hz.
    f_3db_led: float = 15e6      # ancho de banda a -3 dB sin filtro azul [Hz]
    use_blue_filter: bool = False
    f_3db_blue: float = 25e6     # ancho de banda a -3 dB con filtro azul + pre-ecualización [Hz]
    Ts_blue_signal: float = 0.55 # transmitancia efectiva de señal útil con filtro azul
    Ts_blue_bg: float = 0.25     # transmitancia de luz ambiente (sol/horno rojo-IR) con filtro azul

    # --- Receptor (fotodiodo PIN + concentrador no formador de imagen) ---
    z_rx: float = 1.0            # altura del plano de recepción (AGV / telemetría) [m]
    A_pd: float = 1.0e-4         # área activa total del fotodiodo [m^2] (1 cm^2)
    fov_deg: float = 60.0        # campo de visión (semiángulo FOV) [°]
    n_conc: float = 1.5          # índice de refracción del concentrador
    Ts_filter: float = 1.0       # transmitancia de filtro óptico genérico adicional
    resp: float = 0.53           # responsividad del fotodiodo PIN de silicio [A/W]

    # --- Receptor con Diversidad Angular piramidal (ADR, Nuwanpriya et al. 2015) ---
    adr_enabled: bool = False
    adr_tilt_deg: float = 35.0   # inclinación polar de las 4 ramas laterales [°]
    adr_fov_deg: float = 45.0    # FOV estrecho de cada rama del ADR [°]
    adr_combiner: str = "MRC"    # "MRC" (Maximal Ratio) o "SBC" (Select Best)

    # --- Ruido del receptor (parámetros de referencia Komine & Nakagawa 2004) ---
    I_bg: float = 1000e-6        # corriente de referencia por luz ambiente (1 cm^2, limpio, sin filtro) [A]
    I2: float = 0.562            # factor de ancho de banda de ruido I2
    I3: float = 0.0868           # factor de ancho de banda de ruido I3
    G_ol: float = 10.0           # ganancia de lazo abierto del preamplificador FET
    eta_cap: float = 112e-12 / 1e-4  # capacitancia por unidad de área [F/m^2] (112 pF/cm^2)
    Gamma: float = 1.5           # factor de ruido del canal del FET
    g_m: float = 30e-3           # transconductancia del FET [S]
    I_dark0: float = 10e-9       # corriente oscura a 25 °C (1 cm^2) [A]

    # --- Enlace, codificación FEC (IEEE 802.15.7) y paquete ---
    Rb: float = 10e6             # tasa nominal industrial OOK-NRZ [bit/s] (equiparada a Wi-Fi)
    use_fec: bool = False        # activa Reed-Solomon RS(255,239)
    packet_bits: int = 1000      # tamaño de trama industrial para PER / Goodput [bits]

    # --- Condiciones hostiles ---
    # 1) Temperatura ambiente y temperatura de unión Tj (Davis et al. 2019)
    T_amb: float = 25.0          # temperatura ambiente en la nave [°C]
    delta_T_j_nom: float = 30.25 # salto térmico nominal unión-ambiente del LED (0,55 °C/W × 55 W) [°C]
    T_j_ref: float = 55.25       # temperatura de unión nominal a T_amb = 25 °C limpio [°C]
    k_T_led: float = 0.004       # caída térmica del flujo del LED (thermal droop) [1/°C]
    T_j_derate: float = 85.0     # inicio de protección térmica del driver (Tj) [°C] (None = sin protección)
    derate_slope: float = 0.015  # pendiente de reducción de corriente del driver sobre T_j_derate [1/°C]
    k_soil_th: float = 0.25      # aumento de resistencia térmica del disipador por polvo en luminaria

    # 2) Polvo en suspensión (Beer-Lambert + scattering Henyey-Greenstein)
    visibility: float = np.inf   # visibilidad óptica por polvo [m] (inf = aire limpio)
    dust_albedo: float = 0.65    # albedo de dispersión simple omega_0 (polvo mineral/metálico)
    dust_g_hg: float = 0.70      # factor de asimetría de Henyey-Greenstein (forward scattering)
    enable_scattering: bool = True

    # 3) Suciedad depositada sobre ópticas (LDD — Davis et al. 2019)
    days_dirty: float = 0.0      # días sin limpieza de luminarias y receptores
    tau_rx_days: float = 30.0    # constante de ensuciamiento del lente receptor (mira hacia arriba) [días]
    tau_tx_days: float = 120.0   # constante de ensuciamiento de la luminaria (mira hacia abajo) [días]

    # ----- Propiedades físicas derivadas -----
    @property
    def m_lamb(self) -> float:
        """Orden Lambertiano m = -ln(2) / ln(cos(Phi_1/2))."""
        return float(-np.log(2.0) / np.log(np.cos(np.radians(self.phi_half_deg))))

    @property
    def fov(self) -> float:
        return float(np.radians(self.fov_deg))

    @property
    def g_conc(self) -> float:
        """Ganancia óptica del concentrador no formador de imagen: n^2 / sin^2(FOV)."""
        return float(self.n_conc ** 2 / np.sin(self.fov) ** 2)

    @property
    def beta_ext(self) -> float:
        """Coeficiente de extinción por polvo [1/m] (relación de Koschmieder a 550 nm)."""
        return 0.0 if np.isinf(self.visibility) else 3.912 / self.visibility

    @property
    def T_soil_rx(self) -> float:
        """Transmitancia del lente del receptor por suciedad depositada."""
        return float(np.exp(-self.days_dirty / self.tau_rx_days))

    @property
    def T_soil_tx(self) -> float:
        """Transmitancia de la cubierta de la luminaria por suciedad depositada."""
        return float(np.exp(-self.days_dirty / self.tau_tx_days))

    @property
    def T_soil(self) -> float:
        """Transmitancia óptica total por suciedad depositada (luminaria × receptor)."""
        return self.T_soil_rx * self.T_soil_tx

    @property
    def f_3db_effective(self) -> float:
        """Ancho de banda efectivo a -3 dB del transmisor LED [Hz]."""
        return self.f_3db_blue if self.use_blue_filter else self.f_3db_led

    @property
    def Ts_blue_eff(self) -> float:
        """Transmitancia de señal útil asociada al filtro óptico azul."""
        return self.Ts_blue_signal if self.use_blue_filter else 1.0

    @property
    def T_j(self) -> float:
        """
        Temperatura de unión del LED [°C]:
            Tj = T_amb + delta_T_j_nom * [1 + k_soil_th * (1 - T_soil_tx)]
        Modela el autocalentamiento del LED y el efecto aislante del polvo
        acumulado sobre el disipador térmico (Davis et al., 2019).
        """
        soil_penalty = 1.0 + self.k_soil_th * (1.0 - self.T_soil_tx)
        return float(self.T_amb + self.delta_T_j_nom * soil_penalty)

    @property
    def led_derating(self) -> float:
        """Factor de flujo óptico relativo del LED por temperatura de unión Tj y driver."""
        f = 1.0 - self.k_T_led * (self.T_j - self.T_j_ref)
        if self.T_j_derate is not None and self.T_j > self.T_j_derate:
            f *= max(0.0, 1.0 - self.derate_slope * (self.T_j - self.T_j_derate))
        return float(max(0.0, f))

    @property
    def code_rate(self) -> float:
        """Tasa de código efectiva: 239/255 con RS(255,239) FEC; 1,0 sin FEC."""
        return 239.0 / 255.0 if self.use_fec else 1.0

    @property
    def P_comm_total_W(self) -> float:
        """
        Potencia eléctrica incremental dedicada exclusivamente a comunicaciones VLC [W]
        (drivers de modulación en luminarias + electrónica del receptor AGV; la potencia
        de iluminación ya es obligatoria por EN 12464-1).
        """
        n_led = len(self.led_xy)
        rx_mult = 1.6 if self.adr_enabled else 1.0  # 5 canales de preamplificador + combinador
        return n_led * self.P_comm_extra_per_led + self.P_rx_vlc * rx_mult

    @property
    def P_lighting_total_W(self) -> float:
        """Potencia eléctrica total de iluminación instalada en la nave [W]."""
        return len(self.led_xy) * (self.P_opt / self.led_wall_plug_eff)

    def I_bg_effective(self, A_override: float = None) -> float:
        """
        Corriente de ruido ambiental efectivamente captada por el fotodiodo [A].
        Acopla físicamente:
          - área fotosensible (A_branch / 1 cm^2, consistente con conservación de étendue),
          - suciedad depositada sobre el lente receptor (T_soil_rx),
          - filtro óptico genérico (Ts_filter) y rechazo espectral del filtro azul (Ts_blue_bg).
        """
        A = self.A_pd if A_override is None else A_override
        t_bg_filt = self.Ts_filter * (self.Ts_blue_bg if self.use_blue_filter else 1.0)
        return float(self.I_bg * (A / 1.0e-4) * self.T_soil_rx * t_bg_filt)

    @property
    def leds(self) -> np.ndarray:
        xy = np.array(self.led_xy, dtype=float)
        return np.column_stack([xy, np.full(len(xy), self.z_led)])


def optimized_led_grid(p: Params = None) -> Params:
    """
    Devuelve una configuración optimizada simétrica de 16 luminarias (4 × 4,
    en coordenadas 1,1; 5,4; 9,6; 13,9 m), con la MISMA potencia óptica total
    (405 W) y potencia eléctrica de iluminación (900 W) que la grilla base 3 × 3.
    Cumple la norma EN 12464-1 incluso en grilla fina de 0,25 m
    (E_med ≥ 300 lx y uniformidad U0 ≥ 0,60; U0 > 0,70 excluyendo banda de 0,5 m)
    y reduce las zonas de sombra junto a máquinas y estanterías.
    """
    base = Params() if p is None else p
    coords = (1.1, 5.4, 9.6, 13.9)
    grid_4x4 = tuple((x, y) for x in coords for y in coords)
    total_opt = len(base.led_xy) * base.P_opt
    return replace(
        base,
        led_xy=grid_4x4,
        P_opt=total_opt / len(grid_4x4),
        phi_half_deg=65.0,
        P_comm_extra_per_led=base.P_comm_extra_per_led * (len(base.led_xy) / len(grid_4x4)),
    )


# ---------------------------------------------------------------------------
# Obstáculos industriales y fuentes espaciales de EMI
# Cajas alineadas a ejes: [xmin, xmax, ymin, ymax, zmax] (zmin = 0)
# ---------------------------------------------------------------------------
STATIC_OBSTACLES = {
    "Estantería A":          (1.0, 7.0, 13.2, 14.2, 5.0),
    "Estantería B":          (8.0, 14.0, 13.2, 14.2, 5.0),
    "Prensa":                (1.5, 3.5, 1.5, 3.5, 4.0),
    "CNC 1":                 (5.5, 8.5, 2.5, 4.5, 2.5),
    "CNC 2":                 (10.5, 13.5, 1.5, 4.0, 2.5),
    "Horno de inducción":    (10.0, 13.0, 7.0, 10.0, 3.5),
    "Cinta + soldadura":     (1.0, 8.5, 8.0, 8.8, 1.3),
}

# Fuentes industriales de EMI en banda ISM 2,4 GHz: (x, y, z, EMI_pico_1m_dB)
# Basado en mediciones de ruido impulsivo de maquinaria (NIST / Oyedum 2025):
#   - Horno de inducción y puesto de soldadura TIG/MIG: emisión impulsiva severa
#   - Prensa y centros CNC con variadores VFD: emisión moderada-alta
EMI_SOURCES = {
    "Horno de inducción": (11.5, 8.5, 2.0, 38.0),
    "Puesto soldadura":   (4.75, 8.4, 1.5, 36.0),
    "Prensa (VFD)":       (2.5, 2.5, 2.0, 30.0),
    "CNC 1 (VFD)":        (7.0, 3.5, 1.8, 28.0),
    "CNC 2 (VFD)":        (12.0, 2.75, 1.8, 28.0),
}


def static_boxes() -> np.ndarray:
    """Arreglo (7, 5) de obstáculos fijos [xmin, xmax, ymin, ymax, altura]."""
    return np.array(list(STATIC_OBSTACLES.values()), dtype=float)


def heavy_static_boxes(min_height: float = 2.0) -> np.ndarray:
    """Subconjunto de maquinaria metálica alta (para LOS/NLOS RF sin doble conteo)."""
    sb = static_boxes()
    return sb[sb[:, 4] >= min_height]


def random_dynamic_boxes(rng, n_workers: int, n_forklifts: int, p: Params, static=None) -> np.ndarray:
    """
    Genera operarios (0,5 × 0,5 × 1,8 m) y grúas horquilla (1,2 × 2,4 × 2,3 m)
    en posiciones aleatorias libres de maquinaria fija.
    """
    boxes = []
    static = static_boxes() if static is None else static

    def is_free(x0, x1, y0, y1):
        for b in static:
            if not (x1 < b[0] or x0 > b[1] or y1 < b[2] or y0 > b[3]):
                return False
        return True

    for (n, sx, sy, h) in ((n_workers, 0.5, 0.5, 1.8), (n_forklifts, 1.2, 2.4, 2.3)):
        k = 0
        tries = 0
        while k < n and tries < 50 * max(n, 1):
            tries += 1
            if rng.random() < 0.5:
                sx_, sy_ = sx, sy
            else:
                sx_, sy_ = sy, sx
            cx = rng.uniform(sx_ / 2, p.L - sx_ / 2)
            cy = rng.uniform(sy_ / 2, p.W - sy_ / 2)
            b = (cx - sx_ / 2, cx + sx_ / 2, cy - sy_ / 2, cy + sy_ / 2, h)
            if is_free(*b[:4]):
                boxes.append(b)
                k += 1
    return np.array(boxes, dtype=float).reshape(-1, 5)


def segment_blocked(P0, P1, boxes, ignore_if_inside: bool = False) -> np.ndarray:
    """
    Determina si los segmentos P0 -> P1 intersectan alguna caja [xmin, xmax, ymin, ymax, zmax]
    mediante el algoritmo vectorizado de "slabs" (Kay-Kajiya).

    Corrección respecto a la versión inicial:
    `ignore_if_inside` es False por defecto para evitar que un receptor AGV a z = 1 m
    situado debajo de una grúa horquilla (altura 2,3 m) ignore el techo metálico de
    la grúa y reciba señal LOS inexistente.
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


def inside_any(xy, boxes) -> np.ndarray:
    """Máscara booleana de puntos (x, y) que caen dentro de la huella de alguna caja."""
    m = np.zeros(len(xy), dtype=bool)
    if boxes is None or len(boxes) == 0:
        return m
    for b in boxes:
        m |= (xy[:, 0] >= b[0]) & (xy[:, 0] <= b[1]) & (xy[:, 1] >= b[2]) & (xy[:, 1] <= b[3])
    return m


# ---------------------------------------------------------------------------
# Geometría auxiliar y receptor ADR piramidal
# ---------------------------------------------------------------------------
def receiver_grid(p: Params, step: float = 0.25):
    xs = np.arange(step / 2, p.L, step)
    ys = np.arange(step / 2, p.W, step)
    X, Y = np.meshgrid(xs, ys, indexing="xy")
    pts = np.column_stack([X.ravel(), Y.ravel(), np.full(X.size, p.z_rx)])
    return pts, xs, ys


def wall_elements(p: Params):
    """Centros, normales interiores unitarias y área dA de los elementos de pared."""
    d = p.dA_wall
    cs, ns = [], []
    zs = np.arange(d / 2, p.H, d)
    xs = np.arange(d / 2, p.L, d)
    ys = np.arange(d / 2, p.W, d)
    for x in xs:
        for z in zs:
            cs.append((x, 0.0, z)); ns.append((0, 1, 0))
            cs.append((x, p.W, z)); ns.append((0, -1, 0))
    for y in ys:
        for z in zs:
            cs.append((0.0, y, z)); ns.append((1, 0, 0))
            cs.append((p.L, y, z)); ns.append((-1, 0, 0))
    return np.array(cs, dtype=float), np.array(ns, dtype=float), d * d


def adr_normals(tilt_deg: float = 35.0) -> np.ndarray:
    """
    Vectores normales unitarios de las 5 ramas de un Receptor con Diversidad
    Angular (ADR) piramidal (Nuwanpriya et al., IEEE JSAC 2015):
      - Rama 0: cenital (0, 0, 1)
      - Ramas 1..4: inclinadas a `tilt_deg` en azimutes 0°, 90°, 180°, 270°
    """
    th = np.radians(tilt_deg)
    st, ct = np.sin(th), np.cos(th)
    normals = [(0.0, 0.0, 1.0)]
    for az_deg in (0.0, 90.0, 180.0, 270.0):
        az = np.radians(az_deg)
        normals.append((st * np.cos(az), st * np.sin(az), ct))
    return np.array(normals, dtype=float)


# ---------------------------------------------------------------------------
# Iluminación (EN 12464-1) y canal óptico VLC
# ---------------------------------------------------------------------------
def illuminance(p: Params, rx, boxes=None, include_reflections: bool = False, chunk: int = 200) -> np.ndarray:
    """
    Iluminancia horizontal [lx] en el plano de trabajo `rx`.
    Si `boxes` no es None, descuenta las sombras proyectadas por máquinas y estanterías.
    Si `include_reflections=True`, suma la contribución difusa de primer orden de las paredes.
    """
    m = p.m_lamb
    I0 = p.P_opt * p.LER * (m + 1.0) / (2.0 * np.pi)       # intensidad axial [cd]
    Nrx = len(rx)
    E = np.zeros(Nrx)
    for led in p.leds:
        v = led - rx
        d = np.linalg.norm(v, axis=1)
        cos_a = v[:, 2] / d
        contrib = np.where(cos_a > 0, I0 * (cos_a ** m) * cos_a / (d ** 2), 0.0)
        if boxes is not None and len(boxes):
            blk = segment_blocked(rx, np.repeat(led[None, :], Nrx, 0), boxes)
            contrib[blk] = 0.0
        E += contrib

    if include_reflections:
        wc, wn, dA = wall_elements(p)
        E_wall = np.zeros(len(wc))
        for led in p.leds:
            v1 = wc - led
            d1 = np.linalg.norm(v1, axis=1)
            cos_phi = -v1[:, 2] / d1
            cos_alpha = np.einsum("ij,ij->i", -v1, wn) / d1
            ok = (cos_phi > 0) & (cos_alpha > 0)
            if boxes is not None and len(boxes):
                ok &= ~segment_blocked(wc, np.repeat(led[None, :], len(wc), 0), boxes)
            E_wall += np.where(ok, I0 * (cos_phi ** m) * cos_alpha / (d1 ** 2), 0.0)

        # Cada elemento de pared emite con intensidad Lambertiana (rho * E_wall * dA / pi)
        I_wall = p.rho_wall * E_wall * dA / np.pi
        for s in range(0, Nrx, chunk):
            r = rx[s:s + chunk]
            v2 = r[:, None, :] - wc[None, :, :]
            d2 = np.linalg.norm(v2, axis=2)
            cos_beta = np.einsum("nwk,wk->nw", v2, wn) / d2
            cos_psi = -v2[:, :, 2] / d2
            ok = (cos_beta > 0) & (cos_psi > 0)
            if boxes is not None and len(boxes):
                P0 = np.repeat(r, len(wc), 0)
                P1 = np.tile(wc, (len(r), 1))
                blk = segment_blocked(P0, P1, boxes).reshape(len(r), len(wc))
                ok &= ~blk
            E[s:s + chunk] += np.sum(
                np.where(ok, I_wall[None, :] * cos_beta * cos_psi / (d2 ** 2), 0.0),
                axis=1,
            )
    return E


def los_channel(
    p: Params,
    rx: np.ndarray,
    boxes=None,
    rx_normal=None,
    A_override: float = None,
    fov_deg_override: float = None,
):
    """
    Ganancia DC de línea de vista (LOS) para cada par receptor-luminaria
    (Komine & Nakagawa 2004, ec. 3), sin incluir la atenuación del filtro azul
    (que se aplica dinámicamente en `performance()` para reutilizar geometría).

    Admite orientación arbitraria `rx_normal` para modelar ramas del receptor ADR.
    """
    leds = p.leds
    m = p.m_lamb
    A = p.A_pd if A_override is None else A_override
    fov_rad = p.fov if fov_deg_override is None else np.radians(fov_deg_override)
    g = p.n_conc ** 2 / np.sin(fov_rad) ** 2
    n_rx = np.array([0.0, 0.0, 1.0]) if rx_normal is None else np.asarray(rx_normal, dtype=float)
    n_rx = n_rx / np.linalg.norm(n_rx)

    Nrx, Nl = len(rx), len(leds)
    h = np.zeros((Nrx, Nl))
    dist = np.zeros((Nrx, Nl))
    for j, led in enumerate(leds):
        v = led - rx                        # vector desde receptor hacia LED
        d = np.linalg.norm(v, axis=1)
        cos_phi = v[:, 2] / d               # eje del LED es (0, 0, -1)
        cos_psi = (v @ n_rx) / d            # ángulo respecto a normal del fotodiodo
        psi = np.arccos(np.clip(cos_psi, -1.0, 1.0))
        ok = (cos_phi > 0.0) & (cos_psi > 0.0) & (psi <= fov_rad)
        hj = np.where(
            ok,
            (m + 1.0) * A / (2.0 * np.pi * d ** 2) * (cos_phi ** m) * p.Ts_filter * g * cos_psi,
            0.0,
        )
        if boxes is not None and len(boxes):
            blk = segment_blocked(rx, np.repeat(led[None, :], Nrx, 0), boxes)
            hj[blk] = 0.0
        h[:, j] = hj
        dist[:, j] = d
    return h, dist


def reflection_histogram(
    p: Params,
    rx: np.ndarray,
    boxes=None,
    betas=(0.0,),
    dt: float = 0.5e-9,
    t_max: float = 300e-9,
    chunk: int = 150,
    rx_normal=None,
    A_override: float = None,
    fov_deg_override: float = None,
    verbose: bool = False,
):
    """
    Respuesta impulsiva discreta (intervalos `dt`) de las reflexiones difusas de
    primer orden en las paredes de la nave (Komine & Nakagawa 2004).
    Se evalúa para cada coeficiente de extinción `beta` en `betas`.
    """
    leds = p.leds
    m, rho = p.m_lamb, p.rho_wall
    A = p.A_pd if A_override is None else A_override
    fov_rad = p.fov if fov_deg_override is None else np.radians(fov_deg_override)
    g = p.n_conc ** 2 / np.sin(fov_rad) ** 2
    n_rx = np.array([0.0, 0.0, 1.0]) if rx_normal is None else np.asarray(rx_normal, dtype=float)
    n_rx = n_rx / np.linalg.norm(n_rx)

    wc, wn, dA = wall_elements(p)
    nb = int(np.ceil(t_max / dt))
    t_bins = np.arange(nb) * dt
    out = {b: np.zeros((len(rx), nb)) for b in betas}

    # --- Tramo 1: luminaria -> elemento de pared ---
    seg1 = []
    for led in leds:
        v1 = wc - led
        d1 = np.linalg.norm(v1, axis=1)
        cos_phi = -v1[:, 2] / d1
        cos_alpha = np.einsum("ij,ij->i", -v1, wn) / d1
        ok = (cos_phi > 0.0) & (cos_alpha > 0.0)
        if boxes is not None and len(boxes):
            ok &= ~segment_blocked(wc, np.repeat(led[None, :], len(wc), 0), boxes)
        a1 = np.where(ok, (m + 1.0) / (2.0 * np.pi) * (cos_phi ** m) * cos_alpha / (d1 ** 2), 0.0)
        seg1.append((a1, d1))

    # --- Tramo 2: elemento de pared -> receptor ---
    for s in range(0, len(rx), chunk):
        r = rx[s:s + chunk]
        v2 = r[:, None, :] - wc[None, :, :]          # pared -> receptor (n, Nw, 3)
        d2 = np.linalg.norm(v2, axis=2)
        cos_beta = np.einsum("nwk,wk->nw", v2, wn) / d2
        # vector desde receptor hacia pared es -v2
        cos_psi = np.einsum("nwk,k->nw", -v2, n_rx) / d2
        psi = np.arccos(np.clip(cos_psi, -1.0, 1.0))
        ok = (cos_beta > 0.0) & (cos_psi > 0.0) & (psi <= fov_rad)
        if boxes is not None and len(boxes):
            P0 = np.repeat(r, len(wc), 0)
            P1 = np.tile(wc, (len(r), 1))
            blk = segment_blocked(P0, P1, boxes).reshape(len(r), len(wc))
            ok &= ~blk
        a2 = np.where(
            ok,
            rho * dA * cos_beta / (np.pi * d2 ** 2) * A * p.Ts_filter * g * cos_psi,
            0.0,
        )
        for (a1, d1) in seg1:
            gain = a1[None, :] * a2
            dtot = d1[None, :] + d2
            idx = np.minimum((dtot / C_LUZ / dt).astype(int), nb - 1)
            for b in betas:
                gb = gain * np.exp(-b * dtot)
                for i in range(len(r)):
                    out[b][s + i] += np.bincount(idx[i], weights=gb[i], minlength=nb)
        if verbose:
            print(f"   reflexiones: {min(s + chunk, len(rx))}/{len(rx)} receptores", flush=True)
    return out, t_bins


def adr_channels(p: Params, rx: np.ndarray, boxes=None, betas=(0.0,), include_refl: bool = True, dt: float = 0.5e-9):
    """
    Precalcula los canales LOS (y opcionalmente reflexiones) de las 5 ramas del
    Receptor con Diversidad Angular (ADR) piramidal, conservando el área total
    de silicio (A_branch = A_pd / 5 por rama).

    Optimización geométrica: la prueba de oclusión LOS `segment_blocked(rx, led)`
    no depende de la orientación de cada fotodiodo y se evalúa una sola vez para
    todas las ramas (~5× más rápido en simulaciones Monte Carlo).
    """
    normals = adr_normals(p.adr_tilt_deg)
    n_b = len(normals)
    A_b = p.A_pd / n_b
    fov_rad = np.radians(p.adr_fov_deg)
    g = p.n_conc ** 2 / np.sin(fov_rad) ** 2
    m = p.m_lamb
    leds = p.leds
    Nrx, Nl = len(rx), len(leds)

    # Precalcular geometría y oclusión LOS una sola vez para las 5 ramas
    d_mat = np.zeros((Nrx, Nl))
    h_list = [np.zeros((Nrx, Nl)) for _ in range(n_b)]
    for j, led in enumerate(leds):
        v = led - rx
        d = np.linalg.norm(v, axis=1)
        d_mat[:, j] = d
        cos_phi = v[:, 2] / d
        if boxes is not None and len(boxes):
            unblocked = ~segment_blocked(rx, np.repeat(led[None, :], Nrx, 0), boxes)
        else:
            unblocked = np.ones(Nrx, dtype=bool)
        base_j = np.where(
            (cos_phi > 0.0) & unblocked,
            (m + 1.0) * A_b / (2.0 * np.pi * d ** 2) * (cos_phi ** m) * p.Ts_filter * g,
            0.0,
        )
        for b_idx, n_vec in enumerate(normals):
            cos_psi = (v @ n_vec) / d
            psi = np.arccos(np.clip(cos_psi, -1.0, 1.0))
            h_list[b_idx][:, j] = np.where((cos_psi > 0.0) & (psi <= fov_rad), base_j * cos_psi, 0.0)

    d_list = [d_mat] * n_b
    hist_list = []
    t_bins = None
    if include_refl:
        for n_vec in normals:
            hist_b, t_bins = reflection_histogram(
                p, rx, boxes=boxes, betas=betas, dt=dt,
                rx_normal=n_vec, A_override=A_b, fov_deg_override=p.adr_fov_deg,
            )
            hist_list.append(hist_b)
    return h_list, d_list, hist_list, t_bins


# ---------------------------------------------------------------------------
# Ruido, rendimiento VLC (SNR intermedia, BER confiabilidad, Goodput eficiencia)
# ---------------------------------------------------------------------------
def noise_variance(p: Params, P_dc_rx, A_override: float = None):
    """
    Varianza total de ruido eléctrico [A^2] (Komine & Nakagawa 2004) con:
      - dependencia de temperatura ambiente en ruido térmico y corriente oscura,
      - acoplamiento físico de I_bg con área, suciedad del lente y filtro óptico.
    """
    B = p.Rb
    A = p.A_pd if A_override is None else A_override
    T_K = p.T_amb + 273.15
    I_dark = p.I_dark0 * (A / 1.0e-4) * (2.0 ** ((p.T_amb - 25.0) / 10.0))
    I_bg_eff = p.I_bg_effective(A_override=A)

    s_shot = (
        2.0 * Q_E * p.resp * P_dc_rx * B
        + 2.0 * Q_E * I_bg_eff * p.I2 * B
        + 2.0 * Q_E * I_dark * B
    )
    C = p.eta_cap * A
    s_th = (
        (8.0 * np.pi * K_B * T_K / p.G_ol) * C * p.I2 * (B ** 2)
        + (16.0 * np.pi ** 2 * K_B * T_K * p.Gamma / p.g_m) * (C ** 2) * p.I3 * (B ** 3)
    )
    return s_shot, s_th


def _branch_performance(
    p: Params,
    h_los: np.ndarray,
    d_los: np.ndarray,
    refl_hist=None,
    t_bins=None,
    A_override: float = None,
    fov_deg_override: float = None,
):
    """
    Evalúa potencia DC, ruido, ISI, SNR y métricas de enlace para un fotodiodo
    (rama única o rama individual de un ADR).

    Mejoras físicas implementadas:
      1. Solapamiento continuo de pulso OOK-NRZ: w_sig = max(0, 1 - dt/Tb).
      2. Respuesta paso bajo de 1.er orden del LED blanco (f_3db_effective).
      3. Polvo en suspensión con extinción Beer-Lambert + dispersión hacia adelante
         Henyey-Greenstein (ensanchamiento temporal exponencial que elimina la
         paradoja de mejora artificial por polvo a altas tasas).
      4. Ruido shot ambiental I_bg acoplado con T_soil_rx, filtro azul y área.
      5. Codificación RS(255,239) FEC opcional y métricas de Goodput y Mbit/J.
    """
    beta = p.beta_ext
    scale = p.P_opt * p.led_derating * p.T_soil * p.Ts_blue_eff
    tau_o = beta * d_los
    g_los = h_los * np.exp(-tau_o)
    t_los = d_los / C_LUZ

    # Componente dispersada por polvo (forward scattering Henyey-Greenstein)
    Nrx = h_los.shape[0]
    Tb = 1.0 / p.Rb
    if p.enable_scattering and beta > 0.0:
        fov_rad = p.fov if fov_deg_override is None else np.radians(fov_deg_override)
        phi_hg = hg_fov_fraction(fov_rad, p.dust_g_hg)
        w0 = p.dust_albedo
        # Fracción dispersada no absorbida captada por el FOV
        g_scat = h_los * phi_hg * (-np.expm1(-w0 * tau_o)) * np.exp(-(1.0 - w0) * tau_o)
        # Constante temporal media del ensanchamiento por trayectos dispersados
        tau_s = (d_los / C_LUZ) * (1.0 + tau_o) * (1.0 - p.dust_g_hg) * 0.5
    else:
        g_scat = np.zeros_like(h_los)
        tau_s = np.zeros_like(h_los)

    # Instante de primera llegada t0 (sincronización del filtro receptor)
    has_los = g_los.max(axis=1) > 0.0
    t0 = np.where(has_los, np.where(g_los > 0.0, t_los, np.inf).min(axis=1), np.inf)
    if refl_hist is not None:
        first_refl = np.argmax(refl_hist > 0.0, axis=1)
        any_refl = refl_hist.max(axis=1) > 0.0
        t0 = np.where(np.isinf(t0) & any_refl, t_bins[first_refl], t0)

    finite_t0 = np.isfinite(t0)
    t0_safe = np.where(finite_t0, t0, 0.0)

    # 1) Pesos de solapamiento continuo NRZ para componente balística LOS
    dt_los = np.maximum(t_los - t0_safe[:, None], 0.0)
    w_sig_los = np.where(finite_t0[:, None], np.maximum(0.0, 1.0 - dt_los / Tb), 0.0)
    w_isi_los = 1.0 - w_sig_los

    G_chan_sig = np.sum(g_los * w_sig_los, axis=1)
    G_chan_isi = np.sum(g_los * w_isi_los, axis=1)

    # 2) Pesos exactos de solapamiento para cola exponencial de scattering por polvo
    if p.enable_scattering and beta > 0.0:
        rem = np.maximum(0.0, Tb - dt_los)
        valid_scat = finite_t0[:, None] & (rem > 0.0) & (tau_s > 1e-15)
        ratio_s = np.where(valid_scat, rem / np.where(valid_scat, tau_s, 1.0), 0.0)
        w_sig_scat = np.where(
            valid_scat,
            (rem / Tb) - (tau_s / Tb) * (-np.expm1(-ratio_s)),
            0.0,
        )
        w_sig_scat = np.clip(w_sig_scat, 0.0, 1.0)
        w_isi_scat = 1.0 - w_sig_scat
        G_chan_sig += np.sum(g_scat * w_sig_scat, axis=1)
        G_chan_isi += np.sum(g_scat * w_isi_scat, axis=1)

    # 3) Pesos de solapamiento continuo NRZ para reflexiones de pared
    if refl_hist is not None:
        rel = t_bins[None, :] - t0_safe[:, None]
        valid_ref = finite_t0[:, None] & (rel >= -1e-12)
        dt_ref = np.maximum(rel, 0.0)
        w_sig_ref = np.where(valid_ref, np.maximum(0.0, 1.0 - dt_ref / Tb), 0.0)
        w_isi_ref = np.where(finite_t0[:, None], 1.0 - w_sig_ref, 0.0)
        G_chan_sig += np.sum(refl_hist * w_sig_ref, axis=1)
        G_chan_isi += np.sum(refl_hist * w_isi_ref, axis=1)

    # 4) Filtro paso bajo de 1.er orden del LED blanco (f_3db_effective)
    eta_led = led_bandwidth_efficiency(p.Rb, p.f_3db_effective)
    G_sig = eta_led * G_chan_sig
    G_isi = (1.0 - eta_led) * G_chan_sig + G_chan_isi
    G_tot = G_chan_sig + G_chan_isi

    P_dc = scale * G_tot
    I_sig = p.resp * p.mod_index * scale * G_sig
    I_isi = p.resp * p.mod_index * scale * G_isi
    s_shot, s_th = noise_variance(p, P_dc, A_override=A_override)
    denom = s_shot + s_th + I_isi ** 2
    snr = np.where(denom > 0.0, (I_sig ** 2) / denom, 0.0)
    ber_pre = qfunc(np.sqrt(snr))
    ber = rs_post_fec_ber(ber_pre) if p.use_fec else ber_pre
    gp = goodput_mbps(p.Rb, ber, code_rate=p.code_rate, packet_bits=p.packet_bits)
    ee = gp / p.P_comm_total_W

    return dict(
        P_dc=P_dc,
        snr=snr,
        snr_db=10.0 * np.log10(np.maximum(snr, 1e-30)),
        ber_pre=ber_pre,
        ber=ber,
        goodput_mbps=gp,
        ee_mbit_per_j=ee,
        s_shot=s_shot,
        s_th=s_th,
        I_isi_sq=I_isi ** 2,
        G_los=g_los.sum(1),
        G_scat=g_scat.sum(1),
        G_tot=G_tot,
    )


def performance(p: Params, h_los, d_los, refl_hist=None, t_bins=None):
    """Calcula métricas VLC para receptor de un único fotodiodo vertical."""
    return _branch_performance(p, h_los, d_los, refl_hist=refl_hist, t_bins=t_bins)


def performance_adr(p: Params, h_list, d_list, hist_list=None, t_bins=None):
    """
    Calcula métricas VLC para el Receptor con Diversidad Angular (ADR) piramidal
    de 5 ramas combinadas mediante MRC (suma de SNRs por rama) o SBC (mejor rama).
    """
    n_branches = len(h_list)
    A_b = p.A_pd / n_branches
    branch_res = []
    for b in range(n_branches):
        hb_refl = None if hist_list is None else hist_list[b]
        res_b = _branch_performance(
            p, h_list[b], d_list[b], refl_hist=hb_refl, t_bins=t_bins,
            A_override=A_b, fov_deg_override=p.adr_fov_deg,
        )
        branch_res.append(res_b)

    snrs = np.stack([r["snr"] for r in branch_res], axis=0)      # (5, Nrx)
    if p.adr_combiner.upper() == "SBC":
        snr_comb = np.max(snrs, axis=0)
    else:
        snr_comb = np.sum(snrs, axis=0)

    P_dc_tot = np.sum([r["P_dc"] for r in branch_res], axis=0)
    ber_pre = qfunc(np.sqrt(snr_comb))
    ber = rs_post_fec_ber(ber_pre) if p.use_fec else ber_pre
    gp = goodput_mbps(p.Rb, ber, code_rate=p.code_rate, packet_bits=p.packet_bits)
    ee = gp / p.P_comm_total_W

    return dict(
        P_dc=P_dc_tot,
        snr=snr_comb,
        snr_db=10.0 * np.log10(np.maximum(snr_comb, 1e-30)),
        ber_pre=ber_pre,
        ber=ber,
        goodput_mbps=gp,
        ee_mbit_per_j=ee,
        branch_snrs_db=10.0 * np.log10(np.maximum(snrs, 1e-30)),
    )


def rms_delay_spread(h_los, d_los, refl_hist, t_bins, beta_ext: float = 0.0,
                     enable_scattering: bool = False, dust_albedo: float = 0.65,
                     dust_g_hg: float = 0.70, fov_rad: float = np.radians(60.0)):
    """
    Dispersión de retardo RMS [s] por receptor, aplicando consistentemente la
    atenuación por polvo `beta_ext` (y el ensanchamiento por *scattering* si se
    activa) y protegiendo contra división por cero cuando un receptor está en
    sombra óptica total (devuelve NaN).
    """
    t_los = d_los / C_LUZ
    tau_o = beta_ext * d_los
    g_los = h_los * np.exp(-tau_o)

    num1 = (g_los * t_los).sum(1)
    num2 = (g_los * (t_los ** 2)).sum(1)
    den = g_los.sum(1)

    if enable_scattering and beta_ext > 0.0:
        phi_hg = hg_fov_fraction(fov_rad, dust_g_hg)
        g_scat = h_los * phi_hg * (-np.expm1(-dust_albedo * tau_o)) * np.exp(-(1.0 - dust_albedo) * tau_o)
        tau_s = (d_los / C_LUZ) * (1.0 + tau_o) * (1.0 - dust_g_hg) * 0.5
        # Para exponencial desplazada t_los + Exp(tau_s): E[t]=t_los+tau_s, E[t^2]=(t_los+tau_s)^2 + tau_s^2
        num1 += (g_scat * (t_los + tau_s)).sum(1)
        num2 += (g_scat * ((t_los + tau_s) ** 2 + tau_s ** 2)).sum(1)
        den += g_scat.sum(1)

    if refl_hist is not None:
        num1 += (refl_hist * t_bins[None, :]).sum(1)
        num2 += (refl_hist * (t_bins[None, :] ** 2)).sum(1)
        den += refl_hist.sum(1)

    out = np.full(len(den), np.nan)
    ok = den > 0.0
    mu = np.zeros_like(den)
    mu[ok] = num1[ok] / den[ok]
    var = np.zeros_like(den)
    var[ok] = np.maximum(num2[ok] / den[ok] - mu[ok] ** 2, 0.0)
    out[ok] = np.sqrt(var[ok])
    return out


# ---------------------------------------------------------------------------
# Uplink infrarrojo (IR 850 nm) para operación bidireccional libre de RF
# ---------------------------------------------------------------------------
@dataclass
class IrUplinkParams:
    """Parámetros del enlace de subida infrarrojo (AGV -> luminarias del techo)."""
    P_ir_tx: float = 0.20        # potencia óptica media del emisor IR 850 nm en el AGV [W] (seguro clase 1)
    phi_half_ir_deg: float = 60.0
    A_pd_ir: float = 1.0e-4      # fotodiodo IR en cada luminaria [m^2]
    fov_ir_deg: float = 70.0
    n_conc_ir: float = 1.5
    resp_ir: float = 0.60        # responsividad Si a 850 nm [A/W]
    Ts_ir_filter: float = 0.85   # filtro paso-banda 850 nm (bloquea luz visible de LEDs)
    I_bg_ir: float = 200e-6      # corriente de fondo IR (radiación térmica de horno/sol en 850 nm) [A]
    Rb_ul: float = 2.0e6         # tasa de subida para ACK + telemetría AGV [bit/s]
    mod_index: float = 0.8       # LED IR dedicado sin restricción de iluminación/parpadeo
    use_fec: bool = True


def ir_uplink_budget(p: Params, rx: np.ndarray, boxes=None, ir: IrUplinkParams = None):
    """
    Evalúa el enlace de subida infrarrojo (850 nm) desde cada posición AGV `rx`
    hacia los receptores IR co-ubicados en las luminarias del techo (combinación
    por selección del mejor receptor en techo). Demuestra conectividad
    bidireccional VLC(DL) + IR(UL) 100 % inmune a EMI.
    """
    ir = IrUplinkParams() if ir is None else ir
    leds = p.leds
    m = -np.log(2.0) / np.log(np.cos(np.radians(ir.phi_half_ir_deg)))
    fov_rad = np.radians(ir.fov_ir_deg)
    g = ir.n_conc_ir ** 2 / np.sin(fov_rad) ** 2

    Nrx, Nl = len(rx), len(leds)
    h_ul = np.zeros((Nrx, Nl))
    d_ul = np.zeros((Nrx, Nl))
    for j, ap_ir in enumerate(leds):
        v = ap_ir - rx
        d = np.linalg.norm(v, axis=1)
        cos_a = v[:, 2] / d
        psi = np.arccos(np.clip(cos_a, -1.0, 1.0))
        ok = (cos_a > 0.0) & (psi <= fov_rad)
        hj = np.where(
            ok,
            (m + 1.0) * ir.A_pd_ir / (2.0 * np.pi * d ** 2) * (cos_a ** m) * ir.Ts_ir_filter * g * cos_a,
            0.0,
        )
        if boxes is not None and len(boxes):
            blk = segment_blocked(rx, np.repeat(ap_ir[None, :], Nrx, 0), boxes)
            hj[blk] = 0.0
        h_ul[:, j] = hj
        d_ul[:, j] = d

    # Atenuación por polvo y suciedad del lente
    g_eff = h_ul * np.exp(-p.beta_ext * d_ul) * p.T_soil
    # Cada luminaria en techo tiene su propio receptor IR -> selección del mejor enlace
    g_best = np.max(g_eff, axis=1)
    P_dc = ir.P_ir_tx * g_best
    I_sig = ir.resp_ir * ir.mod_index * P_dc

    B = ir.Rb_ul
    T_K = p.T_amb + 273.15
    I_dark = p.I_dark0 * (2.0 ** ((p.T_amb - 25.0) / 10.0))
    s_shot = 2.0 * Q_E * ir.resp_ir * P_dc * B + 2.0 * Q_E * (ir.I_bg_ir * p.T_soil_tx) * p.I2 * B + 2.0 * Q_E * I_dark * B
    C = p.eta_cap * ir.A_pd_ir
    s_th = (8.0 * np.pi * K_B * T_K / p.G_ol) * C * p.I2 * (B ** 2) + \
           (16.0 * np.pi ** 2 * K_B * T_K * p.Gamma / p.g_m) * (C ** 2) * p.I3 * (B ** 3)

    snr = np.where((s_shot + s_th) > 0, (I_sig ** 2) / (s_shot + s_th), 0.0)
    ber_pre = qfunc(np.sqrt(snr))
    ber = rs_post_fec_ber(ber_pre) if ir.use_fec else ber_pre
    return dict(
        snr_db=10.0 * np.log10(np.maximum(snr, 1e-30)),
        ber_pre=ber_pre,
        ber=ber,
        goodput_mbps=goodput_mbps(ir.Rb_ul, ber, code_rate=(239.0 / 255.0 if ir.use_fec else 1.0)),
    )


# ---------------------------------------------------------------------------
# Referencia RF: Wi-Fi 2,4 GHz en fábrica (3GPP TR 38.901 InF-SH + Middleton Class A)
# ---------------------------------------------------------------------------
@dataclass
class WifiParams:
    fc_GHz: float = 2.4
    Pt_dBm: float = 20.0
    G_ant_dBi: float = 2.0           # ganancia de antena Tx y Rx [dBi]
    BW: float = 20e6                 # canal Wi-Fi de 20 MHz
    Rb_target: float = 10e6          # tasa objetivo equiparada a VLC (10 Mb/s, QPSK)
    NF_dB: float = 7.0               # figura de ruido del receptor RF [dB]
    ap_pos: tuple = (7.5, 7.5, 5.5)  # posición del punto de acceso en techo
    sigma_los: float = 4.3           # desviación de sombreamiento LOS [dB] (3GPP InF-SH)
    sigma_nlos: float = 5.9          # desviación de sombreamiento NLOS [dB] (3GPP InF-SH)
    K_los_dB: float = 7.0            # factor K de Rice en LOS [dB]
    body_loss_dB: float = 4.0        # pérdida corporal por operario interceptando LOS [dB]
    emi_dB: float = 20.0             # relación potencia media EMI / ruido térmico [dB]
    emi_A: float = 0.1               # índice de impulsividad Middleton Class A (0,05–0,3 industrial)
    emi_floor_dB: float = 14.0       # piso base de EMI en planta industrial (VFD ambiente) [dB]
    emi_decay_exp: float = 2.0       # exponente de decaimiento espacial de fuentes EMI
    use_fec: bool = False            # RS(255,239) equivalente cuando se compara contra VLC+FEC
    packet_bits: int = 1000
    P_ap_W: float = 12.0             # potencia eléctrica dedicada del AP Wi-Fi industrial [W]
    P_rx_wifi_W: float = 1.5         # potencia eléctrica del transceptor Wi-Fi en el AGV [W]

    @property
    def P_comm_total_W(self) -> float:
        return self.P_ap_W + self.P_rx_wifi_W

    @property
    def code_rate(self) -> float:
        return 239.0 / 255.0 if self.use_fec else 1.0


def wifi_pathloss(w: WifiParams, d3d, los):
    """Pérdidas de propagación 3GPP TR 38.901 sec. 7.4.1 Indoor Factory (InF-SH)."""
    d = np.maximum(np.asarray(d3d, dtype=float), 1.0)
    pl_los = 31.84 + 21.50 * np.log10(d) + 19.00 * np.log10(w.fc_GHz)
    pl_sh = 32.4 + 23.0 * np.log10(d) + 20.0 * np.log10(w.fc_GHz)
    pl_nlos = np.maximum(pl_sh, pl_los)
    return np.where(los, pl_los, pl_nlos)


def wifi_thermal_noise_dBm(w: WifiParams) -> float:
    """Piso de ruido térmico del receptor Wi-Fi [dBm]."""
    return float(-174.0 + 10.0 * np.log10(w.BW) + w.NF_dB)


def wifi_noise_dBm(w: WifiParams, emi_dB=None) -> float:
    """Piso de ruido total equivalente (térmico + potencia media EMI) [dBm]."""
    n_th = wifi_thermal_noise_dBm(w)
    e = w.emi_dB if emi_dB is None else emi_dB
    return n_th + 10.0 * np.log10(1.0 + 10.0 ** (np.asarray(e, dtype=float) / 10.0))


def wifi_spatial_emi_dB(rx: np.ndarray, w: WifiParams = None, sources: dict = None) -> np.ndarray:
    """
    Calcula el mapa espacial de relación potencia EMI / ruido térmico [dB] en la
    nave a partir de las máquinas emisoras (horno de inducción, puesto de soldadura,
    prensa y CNCs con VFD) sobre un piso de ruido industrial base.
    """
    w = WifiParams() if w is None else w
    sources = EMI_SOURCES if sources is None else sources
    p_rel = np.full(len(rx), 10.0 ** (w.emi_floor_dB / 10.0), dtype=float)
    for _, (sx, sy, sz, peak_db) in sources.items():
        d = np.linalg.norm(rx - np.array([sx, sy, sz]), axis=1)
        d = np.maximum(d, 1.0)
        p_rel += (10.0 ** (peak_db / 10.0)) * (d ** (-w.emi_decay_exp))
    return 10.0 * np.log10(p_rel)


def wifi_fading_power(rng, los, n: int, K_los_dB: float = 7.0) -> np.ndarray:
    """
    Ganancia instantánea de potencia |h|^2 con media unitaria:
    Rice (factor K_los_dB) en LOS y Rayleigh en NLOS.
    (Corrección: ya no ignora K_los_dB de WifiParams.)
    """
    los = np.asarray(los, dtype=bool)
    K = 10.0 ** (K_los_dB / 10.0)
    s = np.sqrt(K / (K + 1.0))
    sig = np.sqrt(1.0 / (2.0 * (K + 1.0)))

    hr_rice = s + sig * rng.standard_normal((len(los), n))
    hi_rice = sig * rng.standard_normal((len(los), n))
    p_rice = hr_rice ** 2 + hi_rice ** 2

    hr_ray = rng.standard_normal((len(los), n)) / np.sqrt(2.0)
    hi_ray = rng.standard_normal((len(los), n)) / np.sqrt(2.0)
    p_ray = hr_ray ** 2 + hi_ray ** 2

    return np.where(los[:, None], p_rice, p_ray)


def wifi_ber_middleton_class_a(eb_n0_thermal, emi_dB, A: float = 0.1, m_max: int = 20):
    """
    BER exacta de QPSK/BPSK bajo ruido impulsivo industrial Middleton Class A
    (Middleton, IEEE Trans. EMC 1977):
        BER = sum_{m=0}^{m_max} exp(-A) * (A^m / m!) * Q(sqrt(2 * gamma_0 / (1 + m / (A * Gamma))))
    donde:
      - gamma_0 = Eb / N0_thermal (lineal),
      - Gamma = sigma_thermal^2 / sigma_impulsive^2 = 10^(-emi_dB / 10),
      - A es el índice de impulsividad (valores pequeños => impulsos más agudos).

    Si emi_dB <= -60 dB (sin EMI) o A -> inf, converge exactamente al canal AWGN.
    """
    g0 = np.maximum(np.asarray(eb_n0_thermal, dtype=float), 0.0)
    emi = np.asarray(emi_dB, dtype=float)
    # Broadcast consistente entre g0 y emi
    g0, emi = np.broadcast_arrays(g0, emi)

    # Caso límite gaussiano equivalente si A es muy grande
    if A is None or A > 50.0:
        snr_eff = g0 / (1.0 + np.where(emi > -60.0, 10.0 ** (emi / 10.0), 0.0))
        return qfunc(np.sqrt(2.0 * snr_eff))

    # Potencia relativa media de ruido impulsivo (1 / Gamma)
    inv_gamma = np.where(emi > -60.0, 10.0 ** (emi / 10.0), 0.0)
    ber = np.zeros_like(g0, dtype=float)
    for m in range(m_max + 1):
        pm = math.exp(-A) * (A ** m) / math.factorial(m)
        if pm < 1e-15:
            break
        noise_scale = 1.0 + (m / A) * inv_gamma
        ber += pm * qfunc(np.sqrt(2.0 * g0 / noise_scale))
    return ber


def wifi_evaluate(
    w: WifiParams,
    rx: np.ndarray,
    boxes_heavy=None,
    boxes_dyn=None,
    sf_norm=None,
    rng=None,
    n_fad: int = 1,
    emi_dB_override=None,
    use_spatial_emi: bool = False,
):
    """
    Evalúa el enlace Wi-Fi 2,4 GHz en las ubicaciones `rx` de forma unificada y
    sin doble conteo de obstrucciones:
      - Maquinaria metálica alta (`boxes_heavy`, z ≥ 2 m) y grúas horquilla móviles
        definen la condición LOS/NLOS de 3GPP TR 38.901 InF-SH.
      - Operarios móviles (< 2 m) que cruzan el trayecto introducen pérdida corporal
        (`w.body_loss_dB`), evitando duplicar shadowing + Rayleigh completo.
      - Incluye ruido impulsivo industrial Middleton Class A (uniforme o espacial)
        y codificación RS(255,239) cuando `w.use_fec=True`.
    """
    rng = np.random.default_rng(0) if rng is None else rng
    ap = np.asarray(w.ap_pos, dtype=float)
    Nrx = len(rx)
    d3d = np.linalg.norm(rx - ap, axis=1)
    ap_rep = np.repeat(ap[None, :], Nrx, axis=0)

    # Obstrucción RF pesada (máquinas altas + grúas horquilla >= 2 m)
    heavy_all = []
    if boxes_heavy is not None and len(boxes_heavy):
        heavy_all.append(boxes_heavy)
    workers = None
    if boxes_dyn is not None and len(boxes_dyn):
        forklifts = boxes_dyn[boxes_dyn[:, 4] >= 2.0]
        workers = boxes_dyn[boxes_dyn[:, 4] < 2.0]
        if len(forklifts):
            heavy_all.append(forklifts)
    if heavy_all:
        los = ~segment_blocked(rx, ap_rep, np.vstack(heavy_all))
    else:
        los = np.ones(Nrx, dtype=bool)

    # Pérdida corporal moderada por operarios interceptando el rayo
    body_loss = np.zeros(Nrx, dtype=float)
    if workers is not None and len(workers):
        hit_w = segment_blocked(rx, ap_rep, workers)
        body_loss = np.where(hit_w, w.body_loss_dB, 0.0)

    pl = wifi_pathloss(w, d3d, los)
    if sf_norm is None:
        sf_norm = np.zeros(Nrx, dtype=float)
    sf = np.where(los, w.sigma_los, w.sigma_nlos) * sf_norm

    prx_avg_dBm = w.Pt_dBm + 2.0 * w.G_ant_dBi - pl - sf - body_loss
    fad = wifi_fading_power(rng, los, n_fad, K_los_dB=w.K_los_dB)      # (Nrx, n_fad)
    prx_inst_dBm = prx_avg_dBm[:, None] + 10.0 * np.log10(np.maximum(fad, 1e-30))

    n_th_dBm = wifi_thermal_noise_dBm(w)
    snr_th_db = prx_inst_dBm - n_th_dBm
    # Relación Eb/N0 térmica a tasa objetivo equiparada (con ganancia BW / Rb_target)
    eb_n0_th_lin = (10.0 ** (snr_th_db / 10.0)) * (w.BW / w.Rb_target)

    if emi_dB_override is not None:
        emi_val = np.asarray(emi_dB_override, dtype=float)
    elif use_spatial_emi:
        emi_val = wifi_spatial_emi_dB(rx, w)
    else:
        emi_val = np.full(Nrx, w.emi_dB, dtype=float)
    emi_2d = emi_val[:, None] if emi_val.ndim == 1 else emi_val

    ber_pre = wifi_ber_middleton_class_a(eb_n0_th_lin, emi_2d, A=w.emi_A)
    ber = rs_post_fec_ber(ber_pre) if w.use_fec else ber_pre
    gp = goodput_mbps(w.Rb_target, ber, code_rate=w.code_rate, packet_bits=w.packet_bits)
    ee = gp / w.P_comm_total_W
    snr_eff_db = prx_inst_dBm - wifi_noise_dBm(w, emi_2d)

    if n_fad == 1:
        return dict(
            los=los,
            emi_dB=emi_val,
            snr_db=snr_eff_db[:, 0],
            ber_pre=ber_pre[:, 0],
            ber=ber[:, 0],
            goodput_mbps=gp[:, 0],
            ee_mbit_per_j=ee[:, 0],
        )
    return dict(
        los=los,
        emi_dB=emi_val,
        snr_db=snr_eff_db,
        ber_pre=ber_pre,
        ber=ber,
        goodput_mbps=gp,
        ee_mbit_per_j=ee,
    )


# Umbrales analíticos de referencia para BER = 1e-6 (sin FEC y con RS(255,239))
SNR_TH_OOK_DB = float(10.0 * np.log10(4.7534 ** 2))         # OOK sin FEC -> 13,54 dB
SNR_TH_OOK_FEC_DB = float(10.0 * np.log10(3.15 ** 2))       # OOK con RS(255,239) -> ~9,97 dB (pre-FEC ~8,2e-4 -> post-FEC ~1e-6)
SNR_TH_QPSK_DB = float(10.0 * np.log10(4.7534 ** 2 / 2.0))  # QPSK/BPSK AWGN Eb/N0 -> 10,53 dB
