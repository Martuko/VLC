# -*- coding: utf-8 -*-
"""
tests/test_vlc_core.py
======================
Suite de pruebas unitarias físicas y de regresión para el simulador VLC industrial.
Verifica que todas las correcciones físicas, matemáticas y de comparación con Wi-Fi
se mantengan intactas y libres de artefactos o paradojas.
"""
from dataclasses import replace
from pathlib import Path
import sys
import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from vlc_core import (
    C_LUZ,
    Params,
    WifiParams,
    adr_channels,
    adr_normals,
    goodput_mbps,
    heavy_static_boxes,
    hg_fov_fraction,
    illuminance,
    inside_any,
    ir_uplink_budget,
    led_bandwidth_efficiency,
    los_channel,
    optimized_led_grid,
    performance,
    performance_adr,
    qfunc,
    receiver_grid,
    reflection_histogram,
    rms_delay_spread,
    rs_post_fec_ber,
    segment_blocked,
    static_boxes,
    wifi_ber_middleton_class_a,
    wifi_evaluate,
    wifi_fading_power,
    wifi_spatial_emi_dB,
)


def test_lambertian_on_axis_gain():
    """Ganancia LOS directamente bajo un único LED coincide con la ecuación de Komine (2004)."""
    p = replace(Params(), led_xy=((7.5, 7.5),), z_led=5.8, z_rx=1.0, phi_half_deg=60.0, fov_deg=60.0)
    rx = np.array([[7.5, 7.5, 1.0]])
    h, d = los_channel(p, rx)
    d_exp = 4.8
    m = 1.0  # phi_1/2 = 60° => m = 1
    g_conc = p.n_conc ** 2 / np.sin(np.radians(60.0)) ** 2
    h_exp = (m + 1.0) * p.A_pd / (2.0 * np.pi * d_exp ** 2) * g_conc * p.Ts_filter
    assert np.isclose(d[0, 0], d_exp, rtol=1e-9)
    assert np.isclose(h[0, 0], h_exp, rtol=1e-9)


def test_en12464_optimized_grid_compliance():
    """La grilla optimizada 4×4 conserva potencia total y cumple EN 12464-1 (E≥300 lx, U0≥0,60)."""
    p0 = Params()
    p_opt = optimized_led_grid(p0)
    assert np.isclose(len(p_opt.led_xy) * p_opt.P_opt, len(p0.led_xy) * p0.P_opt, rtol=1e-9)

    rx, _, _ = receiver_grid(p0, 0.25)
    E_base = illuminance(p0, rx, None, include_reflections=True)
    E_opt = illuminance(p_opt, rx, None, include_reflections=True)
    u0_base = E_base.min() / E_base.mean()
    u0_opt = E_opt.min() / E_opt.mean()

    assert E_opt.mean() >= 300.0
    assert u0_opt >= 0.60
    assert u0_opt > u0_base


def test_dust_monotonic_no_paradox():
    """
    Al empeorar la visibilidad por polvo (200 m -> 5 m), la SNR mediana y la
    cobertura deben ser estrictamente no crecientes (sin la paradoja del modelo anterior).
    """
    p0 = Params()
    rx, _, _ = receiver_grid(p0, 1.0)
    h, d = los_channel(p0, rx)
    V_list = [np.inf, 200.0, 100.0, 50.0, 20.0, 10.0, 5.0]
    betas = tuple(sorted(set(0.0 if np.isinf(v) else 3.912 / v for v in V_list)))
    hist, tb = reflection_histogram(p0, rx, None, betas)

    for Rb, blue in ((2e6, False), (10e6, False), (20e6, True), (50e6, True)):
        snrs = []
        for v in V_list:
            p = replace(p0, visibility=v, Rb=Rb, use_blue_filter=blue)
            r = performance(p, h, d, hist[p.beta_ext], tb)
            snrs.append(float(np.median(r["snr_db"])))
        for prev, curr in zip(snrs, snrs[1:]):
            assert curr <= prev + 1e-6, f"Paradoja de polvo en Rb={Rb/1e6} Mb/s: {snrs}"


def test_led_bandwidth_and_blue_filter_tradeoff():
    """
    El LED fosforescente crudo (3 MHz) colapsa a 20 Mb/s, mientras el filtro azul
    (25 MHz) habilita altas tasas y rechaza luz ambiente intensa.
    """
    assert led_bandwidth_efficiency(1e6, 3e6) > 0.99
    assert led_bandwidth_efficiency(20e6, 3e6) < 0.70
    assert led_bandwidth_efficiency(20e6, 25e6) > 0.99

    p0 = Params()
    rx = np.array([[7.5, 7.5, 1.0]])
    h, d = los_channel(p0, rx)
    r_raw = performance(replace(p0, f_3db_led=3e6, use_blue_filter=False, Rb=20e6), h, d)
    r_blue = performance(replace(p0, use_blue_filter=True, Rb=20e6), h, d)
    assert r_blue["snr_db"][0] > r_raw["snr_db"][0] + 10.0


def test_i_bg_coupled_with_soiling_and_filter():
    """La luz ambiente captada debe atenuarse al ensuciarse el receptor o usar filtro azul."""
    p_clean = Params(I_bg=5100e-6, days_dirty=0.0, use_blue_filter=False)
    p_dirty = replace(p_clean, days_dirty=30.0)
    p_blue = replace(p_clean, use_blue_filter=True)

    assert np.isclose(p_clean.I_bg_effective(), 5100e-6, rtol=1e-9)
    assert p_dirty.I_bg_effective() < p_clean.I_bg_effective()
    assert np.isclose(p_dirty.I_bg_effective(), 5100e-6 * p_dirty.T_soil_rx, rtol=1e-9)
    assert np.isclose(p_blue.I_bg_effective(), 5100e-6 * p_blue.Ts_blue_bg, rtol=1e-9)


def test_junction_temperature_and_thermal_derating():
    """La temperatura de unión Tj supera T_amb, aumenta con suciedad y activa derating."""
    p25 = Params(T_amb=25.0, days_dirty=0.0)
    p70 = replace(p25, T_amb=70.0)
    p70_dirty = replace(p70, days_dirty=60.0)

    assert np.isclose(p25.led_derating, 1.0, atol=1e-3)
    assert p70.T_j > 95.0
    assert p70_dirty.T_j > p70.T_j
    assert p70.led_derating < 0.75


def test_segment_blocked_forklift_roof():
    """Un receptor AGV a z = 1 m debajo de una grúa horquilla (z = 2,3 m) debe quedar bloqueado."""
    forklift = np.array([[5.0, 6.2, 5.0, 7.4, 2.3]])
    rx = np.array([[5.6, 6.2, 1.0]])
    led = np.array([[5.6, 6.2, 5.8]])
    assert segment_blocked(rx, led, forklift)[0]


def test_rms_delay_spread_handles_shadow_and_dust():
    """rms_delay_spread no produce ZeroDivisionError en sombra total y responde al polvo."""
    h_zero = np.zeros((1, 2))
    d_dummy = np.array([[5.0, 6.0]])
    refl_zero = np.zeros((1, 10))
    tb = np.arange(10) * 1e-9
    tau_shadow = rms_delay_spread(h_zero, d_dummy, refl_zero, tb)
    assert np.isnan(tau_shadow[0])

    # Con polvo y scattering, la dispersión temporal crece respecto al mismo LOS sin scattering
    h_los = np.array([[1e-5, 0.0]])
    tau_clean = rms_delay_spread(h_los, d_dummy, refl_zero, tb, beta_ext=0.0, enable_scattering=False)
    tau_dust = rms_delay_spread(h_los, d_dummy, refl_zero, tb, beta_ext=0.2, enable_scattering=True)
    assert tau_clean[0] == 0.0
    assert tau_dust[0] > 0.0


def test_rs_fec_and_goodput():
    """RS(255,239) reduce drásticamente BER en torno a 1e-3 y descuenta su overhead del goodput."""
    ber_pre = np.array([1e-5, 5e-4, 1e-3, 1e-2])
    ber_post = rs_post_fec_ber(ber_pre)
    assert ber_post[0] < 1e-18
    assert ber_post[1] < 1e-8
    assert ber_post[2] < 1e-5
    assert ber_post[3] <= ber_pre[3] + 1e-12

    gp_raw = goodput_mbps(10e6, 0.0, code_rate=1.0)
    gp_fec = goodput_mbps(10e6, 0.0, code_rate=239.0 / 255.0)
    assert np.isclose(gp_raw, 10.0)
    assert np.isclose(gp_fec, 10.0 * 239.0 / 255.0)


def test_adr_recovers_blocked_zenith_link():
    """Cuando un obstáculo bloquea la luminaria cenital y las laterales quedan fuera del FOV vertical, el ADR recupera el enlace."""
    p_base = Params(adr_enabled=False, fov_deg=40.0)
    p_adr = replace(p_base, adr_enabled=True, adr_tilt_deg=35.0, adr_fov_deg=45.0)
    normals = adr_normals(35.0)
    assert normals.shape == (5, 3)
    assert np.allclose(np.linalg.norm(normals, axis=1), 1.0)

    # Obstáculo bloquea la luminaria central (7.5, 7.5); la lateral (2.5, 7.5) está a 41,8° (> 40° FOV vertical)
    overhead_blocker = np.array([[7.0, 8.0, 7.0, 8.0, 3.5]])
    rx_edge = np.array([[6.8, 7.5, 1.0]])
    h1, d1 = los_channel(p_base, rx_edge, overhead_blocker)
    r_1pd = performance(p_base, h1, d1)
    ha, da, _, _ = adr_channels(p_adr, rx_edge, overhead_blocker, include_refl=False)
    r_adr = performance_adr(p_adr, ha, da)
    assert r_1pd["snr_db"][0] < 0.0
    assert r_adr["snr_db"][0] > 30.0
    assert r_adr["ber"][0] < 1e-6


def test_wifi_k_factor_middleton_and_spatial_emi():
    """Verifica que K_los_dB afecta wifi_fading_power, Middleton Class A degrada con EMI y el mapa EMI pica cerca del horno."""
    los = np.ones(100, dtype=bool)
    f_low_k = wifi_fading_power(np.random.default_rng(1), los, 500, K_los_dB=0.0)
    f_high_k = wifi_fading_power(np.random.default_rng(1), los, 500, K_los_dB=20.0)
    # Mayor factor K implica menor varianza de potencia alrededor de 1
    assert np.var(f_high_k) < 0.2 * np.var(f_low_k)

    # Middleton Class A converge a AWGN sin EMI y empeora con EMI
    eb_n0 = np.array([100.0])  # 20 dB térmico
    ber_clean = wifi_ber_middleton_class_a(eb_n0, -100.0, A=0.1)
    ber_awgn = qfunc(np.sqrt(2.0 * eb_n0))
    assert np.isclose(ber_clean[0], ber_awgn[0], rtol=1e-9)

    ber_emi20 = wifi_ber_middleton_class_a(eb_n0, 20.0, A=0.1)
    ber_emi35 = wifi_ber_middleton_class_a(eb_n0, 35.0, A=0.1)
    assert ber_emi35[0] > ber_emi20[0] > ber_clean[0]

    # Mapa espacial EMI: más intenso junto al horno (11.5, 8.5) que en la esquina (1.0, 12.0)
    pts = np.array([[11.5, 6.5, 1.0], [1.0, 12.0, 1.0]])
    emi_pts = wifi_spatial_emi_dB(pts)
    assert emi_pts[0] > emi_pts[1] + 4.0


def test_ir_uplink_provides_reliable_return_channel():
    """El enlace de subida IR 850 nm logra alta confiabilidad libre de RF en el área libre."""
    p = replace(optimized_led_grid(Params()), visibility=20.0, days_dirty=7.0, T_amb=70.0)
    sb = static_boxes()
    rx, _, _ = receiver_grid(p, 1.0)
    work = (~inside_any(rx[:, :2], sb)) & (rx[:, 1] <= 13.2)
    ir_res = ir_uplink_budget(p, rx[work], sb)
    cov_ul = 100.0 * np.mean(ir_res["ber"] <= 1e-6)
    assert cov_ul >= 90.0
    assert np.median(ir_res["goodput_mbps"]) >= 1.8


def test_hg_phase_integral_limits():
    """La integral de Henyey-Greenstein vale 0 en FOV=0 y 1 en FOV=pi."""
    assert np.isclose(hg_fov_fraction(0.0, 0.7), 0.0, atol=1e-12)
    assert np.isclose(hg_fov_fraction(np.pi, 0.7), 1.0, atol=1e-12)
    assert 0.0 < hg_fov_fraction(np.radians(60.0), 0.7) < 1.0
