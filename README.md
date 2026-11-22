# Simulación de un Sistema VLC en una Nave Industrial Metalmecánica de Alta EMI

**Proyecto 2 · Comunicaciones por Luz Visible para Industrias Inteligentes**  
**Integrantes:** *(Completar nombres de integrantes)*  
**Semestre:** 2026

---

## Contexto Industrial: ¿Por qué VLC en esta planta?

Este repositorio simula y evalúa un sistema de **Comunicaciones por Luz Visible (VLC)** desplegado en una **nave de manufactura metalmecánica pesada (`15 × 15 × 6 m`)** equipada con horno de inducción, puesto de soldadura de arco, prensa hidráulica, centros de mecanizado CNC con variadores de frecuencia (VFD), estanterías metálicas de 5 m y tránsito de operarios y grúas horquilla.

El sistema comunica en tiempo real consignas de navegación, telemetría y diagnóstico a **10 Mb/s** hacia una flota de vehículos guiados autónomos (**AGVs**) en el plano de `1,0 m`, complementado con un **enlace de subida infrarrojo (IR 850 nm a 2 Mb/s)** para operación bidireccional **100 % libre de RF**:

* **Inmunidad a la EMI industrial:** El horno de inducción, la soldadura y los VFD generan ruido impulsivo en la banda ISM de 2,4 GHz (modelado con **Middleton Class A**), reduciendo la disponibilidad de Wi-Fi hasta `58,3 %` (`76,1 %` con FEC) en la nave, mientras VLC permanece inalterado (`100 %` con receptor ADR + grilla 4×4).
* **Doble uso de la iluminación (Eficiencia energética):** Al modular luminarias LED *high-bay* que ya cumplen la norma **EN 12464-1** (`E ≥ 300 lx`, `U0 ≥ 0,60`), la potencia incremental exclusiva para comunicaciones es baja, superando a Wi-Fi en eficiencia energética (`0,95–0,96 Mbit/J` frente a `0,69 Mbit/J`).
* **Complementariedad frente a polvo extremo:** En atmósferas muy densas (`V < 7–10 m`), la extinción y dispersión (*scattering*) óptica degradan el enlace VLC mientras RF no se ve afectado, justificando arquitecturas híbridas en zonas donde las emisiones RF estén permitidas.

---

## Métricas Oficiales (Enunciado: Confiabilidad y Eficiencia, máx. 2)

1. **Confiabilidad — `BER` (Tasa de Error de Bit, post-FEC cuando aplica):** umbral objetivo industrial `BER ≤ 10⁻⁶`, resumida como **Cobertura espacial** (% del área libre) y **Disponibilidad temporal** (% de ubicaciones con `BER ≤ 10⁻⁶` en ≥ 99 % de realizaciones Monte Carlo).
2. **Eficiencia — `Goodput efectivo [Mb/s]` y `Eficiencia Energética Incremental [Mbit/J]`:** tasa útil libre de errores descontando *overhead* de codificación `RS(255,239)` y pérdida de paquetes (`PER`, trama de 1000 bits), dividida por la potencia eléctrica incremental dedicada a comunicaciones.
* *Variable física intermedia:* La **SNR eléctrica [dB]** se visualiza en mapas de calor para diagnosticar los márgenes de enlace.

---

## Estructura del Repositorio

| Archivo / Carpeta | Descripción |
| --- | --- |
| [vlc_core.py](vlc_core.py) | Núcleo físico: canal Lambertiano LOS + 1.ª reflexión, ancho de banda del LED blanco y filtro azul, polvo (Beer-Lambert + *scattering* Henyey-Greenstein), suciedad óptica LDD, temperatura de unión $T_j$, ruido *shot*/térmico acoplado a $I_{\text{bg}}$, solapamiento NRZ/ISI, **RS(255,239) FEC**, **Receptor ADR piramidal (5 ramas MRC/SBC)**, **Uplink IR 850 nm** y referencia **Wi-Fi 2,4 GHz (3GPP TR 38.901 InF-SH + EMI Middleton Class A)**. |
| [run_all.py](run_all.py) | Ejecuta los 10 experimentos (`E1–E10`) con semilla fija (`2023`) y guarda 13 figuras PNG, 7 archivos CSV, `resumen_resultados.csv` y `log_ejecucion.txt` en `resultados/`. Tarda ≈ 2 min. |
| [demo_interactiva.py](demo_interactiva.py) | Demostración en vivo para presentación oral (~7 s de arranque): controles de polvo, limpieza, temperatura ($T_j$), luz ambiente, operarios, tasa y conmutadores de **Filtro azul**, **RS-FEC**, **ADR + 4×4** y **mapa Wi-Fi bajo EMI**. |
| [tests/test_vlc_core.py](tests/test_vlc_core.py) | Suite de 13 pruebas unitarias físicas (`pytest`) que garantizan monotonicidad ante polvo, cumplimiento EN 12464-1, modelos de LED/ADR/FEC/Middleton y cero regresiones. |
| [docs/MARCO_TEORICO_Y_MODELOS.md](docs/MARCO_TEORICO_Y_MODELOS.md) | Derivación detallada de ecuaciones, contexto industrial y tabla de trazabilidad 1 a 1 a artículos de IEEE Xplore, ScienceDirect, MDPI, Optica y estándares. |
| [docs/REGISTRO_DE_MEJORAS.md](docs/REGISTRO_DE_MEJORAS.md) | Bitácora técnica completa (Antes vs. Después) documentando cada corrección física, metodológica y de software respecto a la versión inicial. |
| `resultados/` | Figuras PNG, barridos CSV, `resumen_resultados.csv` y `log_ejecucion.txt`. |

---

## Instalación y Ejecución

```bash
python3 -m venv .venv
source .venv/bin/activate        # En Windows: .venv\Scripts\activate
pip install -r requirements.txt

# 1. Ejecutar suite de verificación física (requiere pytest: pip install pytest)
pytest -v tests/test_vlc_core.py

# 2. Regenerar todos los experimentos E1–E10 (≈ 2 min)
python run_all.py

# 3. Abrir demostración interactiva en vivo (entorno gráfico)
python demo_interactiva.py
```

---

## Experimentos y Resultados Principales (`10 Mb/s`)

| Id | Experimento | Figura(s) | Hallazgo clave |
| --- | --- | --- | --- |
| **E1** | Planta EMI e iluminancia EN 12464-1 | `E1_planta_iluminancia.png` | Grilla base `3×3` incumple uniformidad (`349 lx`, `U0 = 0,46`); grilla optimizada `4×4` a igual potencia cumple EN 12464-1 (`302 lx`, `U0 = 0,61`). |
| **E2** | Caso de referencia (`10 Mb/s`) | `E2_referencia_potencia_snr_retardo.png`, `E2_respuesta_impulsiva.png` | Potencia `-8,5` a `-3,6 dBm`, SNR `17,4` a `34,5 dB` (med. `28,7 dB`), Goodput `10,00 Mb/s`, $\tau_{\text{RMS}} \le 7,3\text{ ns}$. |
| **E3** | Seis escenarios hostiles (`S0–S5`) | `E3_mapas_snr_escenarios.png`, `E3_cdf_ber.png` | Factores aislados mantienen `100 %`; obstáculos (bloqueando LOS y reflexiones) bajan a `92,6 %` (`9,46 Mb/s`); combinado `S5` cae a `87,0 %` (`8,85 Mb/s`). |
| **E4** | Polvo en suspensión (extinción + *scattering*) | `E4_polvo_ber.png` | Degradación estrictamente monótona sin paradojas: a `10 Mb/s` soporta hasta `V = 7 m` (`100 %`) y cae a `59,1 %` en `V = 5 m`. |
| **E5** | Suciedad en ópticas (LDD + acople $T_j$/$I_{\text{bg}}$) | `E5_suciedad_limpieza.png` | Con polvo intenso (`V = 10 m`) el ensuciamiento acelera la caída de cobertura tras 60–100 días; limpieza semanal/mensual evita degradación acumulada. |
| **E6** | Temperatura ambiente y unión $T_j$ | `E6_temperatura.png` | Sobre $T_{\text{amb}} \approx 53\,^\circ\text{C}$ ($T_j \approx 85\,^\circ\text{C}$) se activa el *derating* del driver (hasta `-1,0` a `-1,4 dB` en régimen limitado por ruido); la corriente oscura es despreciable. |
| **E7** | Tasa de bits, BW del LED e ISI | `E7_tasa_bits.png` | LED crudo (`3 MHz`) colapsa sobre `5–10 Mb/s`; pre-ecualizado (`15 MHz`) llega a `30 Mb/s`; filtro azul (`25 MHz`) alcanza `41+ Mb/s` ideal y `24,5 Mb/s` en escenario hostil con FEC. |
| **E8** | Obstáculos móviles (Monte Carlo equiparado) | `E8_obstaculos_montecarlo.png` | Con 40 operarios + 4 grúas, VLC 1-PD cae a `70 %` y Wi-Fi bajo EMI espacial a `48 %`; **VLC con ADR piramidal mantiene `~100 %`**. |
| **E9** | Comparación VLC vs Wi-Fi (EMI, polvo, mapa y eficiencia) | `E9_comparacion_wifi.png`, `E9_mapa_espacial_eficiencia.png` | Bajo EMI espacial de planta, VLC ADR+4×4 logra `100,0 %` (`0,96 Mbit/J`) vs Wi-Fi+FEC `76,1 %` (`0,69 Mbit/J`). **Uplink IR 850 nm**: `98,9 %` cobertura (`1,86 Mb/s`). |
| **E10** | Estudio de ablación en escenario combinado hostil | `E10_diseno_mejorado.png` | Disponibilidad ≥ 99 % sube de `60,0 %` (Base) → `62,2 %` (+Filtro+7d) → `71,7 %` (+RS-FEC) → **`78,9 %` (VLC Autónomo ADR+4×4, `0,95 Mbit/J`)**, superando a Wi-Fi+FEC (`70,6 %`, `0,69 Mbit/J`). Modo híbrido alcanza `99,4 %` (`0,40 Mbit/J`). |

---

## Referencias Científicas Principales (Verificadas)

1. **T. Komine y M. Nakagawa**, *"Fundamental analysis for visible-light communication system using LED lights"*, *IEEE Trans. Consumer Electronics*, vol. 50, no. 1, pp. 100–107, 2004. DOI: `10.1109/TCE.2004.1277847` *(IEEE Xplore)*.
2. **M. Uysal, F. Miramirkhani et al.**, *"IEEE 802.15.7r1 reference channel models for visible light communications"*, *IEEE Communications Magazine*, vol. 55, no. 1, pp. 212–217, 2017. DOI: `10.1109/MCOM.2017.1600872CM` *(IEEE Xplore)*.
3. **H. Le Minh et al.**, *"100-Mb/s NRZ Visible Light Communications Using a Postequalized White LED"*, *IEEE Photonics Technology Letters*, vol. 21, no. 15, pp. 1063–1065, 2009. DOI: `10.1109/LPT.2009.2022413` *(IEEE Xplore)*.
4. **A. Nuwanpriya, S.-W. Ho y C. S. Chen**, *"Indoor MIMO Visible Light Communications: Novel Angle Diversity Receivers for Mobile Users"*, *IEEE JSAC*, vol. 33, no. 9, pp. 1780–1792, 2015. DOI: `10.1109/JSAC.2015.2432514` *(IEEE Xplore)*.
5. **A. Al-Kinani et al.**, *"Optical Wireless Communication Channel Measurements and Models"*, *IEEE Commun. Surveys & Tutorials*, vol. 20, no. 3, pp. 1939–1962, 2018. DOI: `10.1109/COMST.2018.2838096` *(IEEE Xplore)*.
6. **L. E. M. Matheus et al.**, *"Visible Light Communication: Concepts, Applications and Challenges"*, *IEEE Commun. Surveys & Tutorials*, vol. 21, no. 4, pp. 3204–3237, 2019. DOI: `10.1109/COMST.2019.2913348` *(IEEE Xplore)*.
7. **P. Palacios Játiva et al.**, *"A VLC Channel Model for Underground Mining Environments With Scattering and Shadowing"*, *IEEE Access*, vol. 8, pp. 185445–185464, 2020 (`10.1109/ACCESS.2020.3030615`), y *"An Enhanced VLC Channel Model for Underground Mining Environments Considering a 3D Dust Particle Distribution Model"*, *MDPI Sensors*, vol. 22, no. 7, 2483, 2022 (`10.3390/s22072483`).
8. **D. Schneider et al.**, *"Impact of industrial environments on visible light communication"*, *Optics Express*, vol. 29, no. 11, 2021. DOI: `10.1364/OE.421757` *(Optica)*.
9. **C. Jenila y R. K. Jeyachitra**, *"Green indoor optical wireless communication systems: Pathway towards pervasive deployment"*, *Digital Communications and Networks*, vol. 7, no. 3, pp. 410–444, 2021. DOI: `10.1016/j.dcan.2020.09.004` *(ScienceDirect)*.
10. **D. Middleton**, *"Statistical-Physical Models of Electromagnetic Interference"*, *IEEE Trans. EMC*, vol. EMC-19, no. 3, pp. 106–127, 1977 (`10.1109/TEMC.1977.303527`), y **J.-D. C. Oyedum et al.**, *"Measurement, Modeling and Analysis of the Impact of Machine-Generated Noise on OFDM Signal Reception in Advanced Manufacturing"*, *IEEE GLOBECOM*, 2025 (`10.1109/GLOBECOM59602.2025.11432514`).
11. **R. G. Davis, A. Wilkerson y B. R. Kinzey**, *"Luminaire dirt depreciation (LDD): Field data from several exterior lighting projects"*, *LEUKOS*, vol. 15, no. 1, pp. 55–63, 2019. DOI: `10.1080/15502724.2018.1517598`.
12. **Estándares:** IEEE Std 802.15.7-2018 (`10.1109/IEEESTD.2019.8697198`), 3GPP TR 38.901 (sec. 7.4.1 InF-SH) y EN 12464-1.
