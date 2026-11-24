# Simulación de un Sistema VLC en una Nave Industrial Metalmecánica de Alta EMI

**Proyecto 2 · Comunicaciones por Luz Visible para Industrias Inteligentes**  
**Integrantes:** *(Completar nombres de integrantes)*  
**Semestre:** 2026

**Simulador 3D en el navegador:** abre `docs/index.html` o publícalo con GitHub Pages (ver más abajo). Permite mover un AGV por la nave, cambiar polvo, calor, luz ambiente, operarios, tasa, filtro, FEC y receptor, desplegar un **panel dinámico de gráficos (`E2–E11`)** debajo del modelo 3D y comparar en vivo con Wi-Fi bajo interferencia electromagnética.

> **Revisión técnica 2 (noviembre 2026).** Se corrigió el modelo de interferencia intersimbólica (ahora por apertura de ojo con muestreo óptimo), que dejaba el enlace saturado por ISI; se reemplazó una referencia no verificable; se agregó el experimento E11 y el simulador 3D web. Detalle completo en [docs/REVISION_TECNICA.md](docs/REVISION_TECNICA.md).

---

## Contexto Industrial: ¿Por qué VLC en esta planta?

Este repositorio simula y evalúa un sistema de **Comunicaciones por Luz Visible (VLC)** desplegado en una **nave de manufactura metalmecánica pesada (`15 × 15 × 6 m`)** equipada con horno de inducción, puesto de soldadura de arco, prensa hidráulica, centros de mecanizado CNC con variadores de frecuencia (VFD), estanterías metálicas de 5 m y tránsito de operarios y grúas horquilla.

El sistema comunica en tiempo real consignas de navegación, telemetría y diagnóstico a **10 Mb/s** hacia una flota de vehículos guiados autónomos (**AGVs**) en el plano de `1,0 m`, complementado con un **enlace de subida infrarrojo (IR 850 nm a 2 Mb/s)** para operación bidireccional **100 % libre de RF**:

* **Inmunidad a la EMI industrial:** El horno de inducción, la soldadura y los VFD generan ruido impulsivo en la banda ISM de 2,4 GHz (modelado con **Middleton Class A**), reduciendo la disponibilidad de Wi-Fi a `58,3 %` (`76,1 %` con FEC) en la nave, mientras VLC no se ve afectado (`99,4 %` con la instalación base de 9 luminarias y `100 %` con ADR + grilla 4×4, en aire limpio).
* **Doble uso de la iluminación (Eficiencia energética):** Al modular luminarias LED *high-bay* (la grilla 4×4 cumple **EN 12464-1** en la nave vacía: `302 lx`, `U0 = 0,61`), la potencia incremental para comunicar es baja. VLC supera a Wi-Fi en Mbit/J **solo si** el driver de modulación consume poco: con un AP de 12 W, menos de `1,56 W` por luminaria en la grilla 3×3 o `0,80 W` en la 4×4 (experimento E11).
* **Riesgo principal de VLC: el bloqueo.** El polvo, el calor y la luz solar consumen margen (5–15 dB) pero no cortan el enlace a 10 Mb/s; lo que sí lo corta es la sombra de estanterías, máquinas y grúas, y la suciedad acumulada sin limpieza. Por eso conviene un respaldo RF donde las emisiones estén permitidas (híbrido: `99,4 %`).

---

## Métricas Oficiales (Enunciado: Confiabilidad y Eficiencia, máx. 2)

1. **Confiabilidad — `BER` (Tasa de Error de Bit, post-FEC cuando aplica):** umbral objetivo industrial `BER ≤ 10⁻⁶`, resumida como **Cobertura espacial** (% del área libre) y **Disponibilidad temporal** (% de ubicaciones con `BER ≤ 10⁻⁶` en ≥ 99 % de realizaciones Monte Carlo).
2. **Eficiencia — `Goodput efectivo [Mb/s]` y `Eficiencia Energética Incremental [Mbit/J]`:** tasa útil libre de errores descontando *overhead* de codificación `RS(255,239)` y pérdida de paquetes (`PER`, trama de 1000 bits), dividida por la potencia eléctrica incremental dedicada a comunicaciones.
* *Variable física intermedia:* la **SNR de apertura de ojo [dB]** (peor caso de ISI) se visualiza en mapas de calor para diagnosticar los márgenes de enlace.

---

## Estructura del Repositorio

| Archivo / Carpeta | Descripción |
| --- | --- |
| [vlc_core.py](vlc_core.py) | Núcleo físico: canal Lambertiano LOS + 1.ª reflexión, ancho de banda del LED blanco y filtro azul, polvo (Beer-Lambert + *scattering* Henyey-Greenstein), suciedad óptica LDD, temperatura de unión $T_j$, ruido *shot*/térmico acoplado a $I_{\text{bg}}$, solapamiento NRZ/ISI, **RS(255,239) FEC**, **Receptor ADR piramidal (5 ramas MRC/SBC)**, **Uplink IR 850 nm** y referencia **Wi-Fi 2,4 GHz (3GPP TR 38.901 InF-SH + EMI Middleton Class A)**. |
| [run_all.py](run_all.py) | Ejecuta los 11 experimentos (`E1–E11`) con semilla fija (`2023`) y guarda 14 figuras PNG, 8 archivos CSV (incluido `resumen_resultados.csv`) y `log_ejecucion.txt` en `resultados/`. Tarda ≈ 3 min. |
| [demo_interactiva.py](demo_interactiva.py) | Demostración en vivo para presentación oral (~7 s de arranque): controles de polvo, limpieza, temperatura ($T_j$), luz ambiente, operarios, tasa y conmutadores de **Filtro azul**, **RS-FEC**, **ADR + 4×4** y **mapa Wi-Fi bajo EMI**. |
| [tests/test_vlc_core.py](tests/test_vlc_core.py) | Suite de 16 pruebas unitarias físicas (`pytest`): monotonicidad ante polvo, EN 12464-1, LED/ADR/FEC/Middleton, y (revisión 2) que los obstáculos nunca mejoren la SNR y que calor y luz ambiente sí la reduzcan. |
| [docs/index.html](docs/index.html) | **Simulador 3D web** autocontenido (Three.js + panel desplegable de gráficos dinámicos `E2–E11` en Canvas 2D). Se construye con `python web/build.py` desde `web/`. |
| [web/](web/) | Fuente del simulador web: `vlc_physics.js` (port de `vlc_core.py`), `vlc_charts.js` (gráficos dinámicos `E2–E11` sin dependencias), `index.template.html`, `build.py` y las validaciones contra Python (`generar_referencia.py`, `validar_contra_python.js`, `validar_graficos.js`). |
| [docs/REVISION_TECNICA.md](docs/REVISION_TECNICA.md) | Revisión técnica 2: hallazgos, correcciones, verificación de referencias y resultados antes/después. |
| [docs/MARCO_TEORICO_Y_MODELOS.md](docs/MARCO_TEORICO_Y_MODELOS.md) | Derivación detallada de ecuaciones, contexto industrial y tabla de trazabilidad 1 a 1 a artículos de IEEE Xplore, ScienceDirect, MDPI, Optica y estándares. |
| [docs/REGISTRO_DE_MEJORAS.md](docs/REGISTRO_DE_MEJORAS.md) | Bitácora técnica completa (Antes vs. Después) documentando cada corrección física, metodológica y de software respecto a la versión inicial. |
| `resultados/` | Figuras PNG, barridos CSV, `resumen_resultados.csv` y `log_ejecucion.txt`. |

---

## Instalación y Ejecución

```bash
python3 -m venv .venv
source .venv/bin/activate        # En Windows: .venv\Scripts\activate
source .venv/bin/activate.fish   # En terminal fish

pip install -r requirements.txt

# 1. Ejecutar suite de verificación física (requiere pytest: pip install pytest)
pytest -v tests/test_vlc_core.py

# 2. Regenerar todos los experimentos E1–E11 (≈ 3 min)
python run_all.py

# 3. Abrir demostración interactiva en vivo (entorno gráfico)
python demo_interactiva.py

# 4. Simulador 3D web: reconstruir y validar contra Python (requiere Node.js)
python web/build.py
python web/generar_referencia.py && node web/validar_contra_python.js && node web/validar_graficos.js
```

### Publicar el simulador 3D con GitHub Pages

1. Sube el repositorio con la carpeta `docs/` (incluye `index.html` y `.nojekyll`).
2. En GitHub: **Settings → Pages → Build and deployment → Source: Deploy from a branch**.
3. Elige la rama `main` y la carpeta **`/docs`**, y guarda.
4. Tras uno o dos minutos la página queda en `https://TU_USUARIO.github.io/NOMBRE_DEL_REPO/`.

La página no necesita servidor: también funciona abriendo `docs/index.html` directamente en el navegador (con conexión a internet para cargar Three.js).

---

## Experimentos y Resultados Principales (`10 Mb/s`, modelo de ojo)

| Id | Experimento | Figura(s) | Hallazgo clave |
| --- | --- | --- | --- |
| **E1** | Planta EMI e iluminancia EN 12464-1 | `E1_planta_iluminancia.png` | Grilla base `3×3`: `349 lx`, `U0 = 0,46` (no cumple); grilla `4×4` a igual potencia: `302 lx`, `U0 = 0,61` en la nave vacía. Con máquinas, la uniformidad en la zona de tránsito es `0,08` (3×3) y `0,05` (4×4). |
| **E2** | Caso de referencia | `E2_referencia_potencia_snr_retardo.png`, `E2_respuesta_impulsiva.png` | Potencia `-8,5` a `-3,6 dBm`; SNR `49,5` a `58,3 dB` (mediana `56,2 dB`); goodput `10,00 Mb/s`; $\tau_{\text{RMS}} \le 7,3\text{ ns}$. |
| **E3** | Seis escenarios hostiles (`S0–S5`) | `E3_mapas_snr_escenarios.png`, `E3_cdf_ber.png` | Cada factor aislado mantiene `100 %` pero consume margen (polvo + suciedad `-11,3 dB`, sol `-6,0 dB`, calor `-3,7 dB`); obstáculos `95,8 %`; combinado `S5` `89,7 %` (`9,11 Mb/s`). |
| **E4** | Polvo en suspensión (extinción + *scattering*) | `E4_polvo_ber.png` | Degradación monótona; a `10 Mb/s` la cobertura se mantiene en `100 %` incluso con `V = 5 m` (SNR mediana `40,7 dB`). |
| **E5** | Suciedad en ópticas | `E5_suciedad_limpieza.png` | Con polvo denso (`V = 10 m`) la cobertura cae de `90,8 %` a `55,2 %` a los 120 días sin limpieza; en aire limpio, a `85,1 %`. |
| **E6** | Temperatura ambiente y de unión $T_j$ | `E6_temperatura.png` | `-4,3 dB` a `70 °C` y `-8,0 dB` a `85 °C` (incluye el *derating* del driver sobre $T_j = 85\,^\circ\text{C}$); la corriente oscura es despreciable. |
| **E7** | Tasa de bits, BW del LED e ISI | `E7_tasa_bits.png` | Tasa máxima con `100 %` de cobertura (nave ideal): LED sin ecualizar (`3 MHz`) `20 Mb/s`, pre-ecualizado (`15 MHz`) `50 Mb/s`, filtro azul (`25 MHz`) `75 Mb/s`. |
| **E8** | Obstáculos móviles (Monte Carlo) | `E8_obstaculos_montecarlo.png` | Con 40 operarios: fotodiodo FOV 60° `95,6 %`, **fotodiodo FOV 80° `100 %`**, 4 luminarias `84,4 %`, ADR `100 %`, Wi-Fi bajo EMI `48,3 %`. |
| **E9** | VLC vs Wi-Fi (EMI, polvo, mapa, eficiencia) | `E9_comparacion_wifi.png`, `E9_mapa_espacial_eficiencia.png` | Bajo EMI espacial: VLC base `99,4 %`, VLC ADR+4×4 `100 %`, Wi-Fi `58,3 %` (`76,1 %` con FEC). Uplink IR 850 nm: `98,9 %` (`1,86 Mb/s`). |
| **E10** | Ablación en escenario combinado hostil | `E10_diseno_mejorado.png` | Disponibilidad ≥ 99 %: Base `62,8 %` → +Filtro +7 d `69,4 %` → +RS-FEC `79,4 %`; ADR+4×4 `78,9 %` (no mejora sobre FEC); Wi-Fi+FEC `70,6 %`; híbrido `99,4 %`. |
| **E11** | Sensibilidad de la eficiencia (revisión 2) | `E11_sensibilidad_eficiencia.png` | VLC es más eficiente que un AP de 12 W solo si el driver consume menos de `1,56 W` (3×3) o `0,80 W` (4×4) por luminaria; con un AP de 6 W, `0,84 W` y `0,42 W`. |

**Conclusión.** VLC entrega 10 Mb/s confiables en casi toda la nave aunque haya polvo, calor y luz solar, y no se ve afectado por la EMI que degrada a Wi-Fi. Sus límites son el bloqueo de línea de vista y la suciedad acumulada, que se mitigan con más luminarias, un FOV amplio, limpieza periódica y un respaldo RF.

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
10. **D. Middleton**, *"Statistical-Physical Models of Electromagnetic Interference"*, *IEEE Trans. EMC*, vol. EMC-19, no. 3, pp. 106–127, 1977 (`10.1109/TEMC.1977.303527`), y **K. L. Blackard, T. S. Rappaport y C. W. Bostian**, *"Measurements and models of radio frequency impulsive noise for indoor wireless communications"*, *IEEE JSAC*, vol. 11, no. 7, pp. 991–1001, 1993 (`10.1109/49.233212`). Los valores pico de las fuentes EMI de la nave son supuestos de diseño.
11. **R. G. Davis, A. Wilkerson y B. R. Kinzey**, *"Luminaire dirt depreciation (LDD): Field data from several exterior lighting projects"*, *LEUKOS*, vol. 15, no. 1, pp. 55–63, 2019. DOI: `10.1080/15502724.2018.1517598`. Se cita por el concepto de LDD; sus mediciones son de luminarias de exterior, y las constantes de ensuciamiento y el acoplamiento térmico del simulador son supuestos para polvo industrial denso.
12. **Estándares:** IEEE Std 802.15.7-2018 (`10.1109/IEEESTD.2019.8697198`), 3GPP TR 38.901 (sec. 7.4.1 InF-SH) y EN 12464-1.
