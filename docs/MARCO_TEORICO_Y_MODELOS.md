# Marco Teórico, Contexto Industrial y Modelos Físicos

Este documento detalla el contexto de aplicación industrial, el fundamento matemático-físico de cada módulo implementado en [vlc_core.py](../vlc_core.py) y la trazabilidad exacta de todos los parámetros hacia artículos revisados por pares de **IEEE Xplore**, **ScienceDirect**, **MDPI**, **Optica** y estándares internacionales (**IEEE 802.15.7**, **3GPP TR 38.901**, **EN 12464-1**).

---

## 1. Contexto Industrial y Motivación Técnica de VLC

### 1.1 Escenario de aplicación: Nave de manufactura metalmecánica pesada con alta EMI
El simulador modela una nave de manufactura pesada de **15 × 15 × 6 m** con siete equipos/estructuras fijas y tránsito continuo de personal y maquinaria móvil:

| Elemento industrial | Posición `(x, y)` [m] | Altura [m] | Rol físico en el simulador |
| --- | --- | --- | --- |
| **Estanterías metálicas A y B** | `[1–7, 13,2–14,2]` y `[8–14, 13,2–14,2]` | 5,0 | Bloqueo severo de línea de vista (LOS) y reflexiones cerca de la pared norte. |
| **Prensa hidráulica/eléctrica** | `[1,5–3,5, 1,5–3,5]` | 4,0 | Obstáculo alto y fuente de EMI por variadores de frecuencia (VFD, pico 30 dB a 1 m). |
| **Centros CNC 1 y CNC 2** | `[5,5–8,5, 2,5–4,5]` y `[10,5–13,5, 1,5–4,0]` | 2,5 | Obstáculos metálicos y fuentes de ruido impulsivo por VFD (pico 28 dB a 1 m). |
| **Horno de inducción** | `[10,0–13,0, 7,0–10,0]` | 3,5 | Fuente de calor extremo ($T_{\text{amb}}$ hasta 70 °C), radiación óptica de fondo y EMI severa (pico 38 dB a 1 m). |
| **Cinta transportadora + soldadura** | `[1,0–8,5, 8,0–8,8]` | 1,3 | Obstáculo bajo y puesto de soldadura de arco TIG/MIG (EMI impulsiva severa, pico 36 dB a 1 m). |
| **Operarios y grúas horquilla** | Dinámicos (Monte Carlo) | 1,8 / 2,3 | Bloqueos aleatorios transitorios en la zona de tránsito de AGVs ($z_{\text{rx}} = 1,0\text{ m}$). |

### 1.2 ¿Por qué VLC es más conveniente que Wi-Fi en este entorno?
1. **Inmunidad total a la Interferencia Electromagnética (EMI):** En plantas con hornos de inducción, arcos de soldadura TIG/MIG y motores con VFD, las conmutaciones de alta potencia generan ruido impulsivo de banda ancha que se derrama sobre la banda ISM de 2,4 GHz (Middleton, 1977; Blackard et al., 1993). Esto degrada drásticamente la tasa de error y disponibilidad de Wi-Fi, mientras que el enlace óptico VLC opera a longitudes de onda de 400–700 nm sin acoplamiento electromagnético alguno.
2. **Doble uso de la infraestructura de iluminación (Eficiencia energética):** La norma europea **EN 12464-1** exige una iluminancia mínima de $\ge 300\text{ lx}$ y uniformidad $U_0 \ge 0,60$ en áreas de producción industrial. Al modular las luminarias LED *high-bay* ya necesarias para iluminar la planta, el consumo energético incremental dedicado exclusivamente a comunicaciones se reduce a la electrónica de modulación y al receptor óptico (~10–12 W totales para toda la nave), pudiendo lograr mayor eficiencia energética de comunicación (`Mbit/J`) que un punto de acceso RF dedicado siempre que el driver de modulación consuma poco (ver E11: con un AP de 12 W, menos de 1,56 W por luminaria en la grilla 3×3 y 0,80 W en la 4×4), y cuando el Wi-Fi sufre retransmisiones por EMI (Matheus et al., 2019; Jenila & Jeyachitra, 2021).
3. **Confinamiento físico y seguridad de capa 1:** La luz visible no atraviesa muros opacos, confinando las señales de control y telemetría de AGVs dentro del recinto fabril sin fugas hacia el exterior ni congestión del espectro RF licenciado/no licenciado.
4. **Arquitectura bidireccional libre de RF (`VLC Downlink + IR 850 nm Uplink`):** Para mantener inmunidad total a EMI en ambos sentidos, el enlace descendente VLC (`10 Mb/s` nominales para consignas de navegación, mapas y diagnóstico hacia AGVs) se complementa con un enlace de subida infrarrojo a `850 nm` (`2 Mb/s` para telemetría de retorno y ACKs), eliminando la dependencia de RF en zonas críticas.
5. **Complementariedad honesta ante polvo extremo:** Cuando la concentración de polvo en suspensión reduce la visibilidad por debajo de ~7–10 m, la extinción y dispersión óptica degradan al VLC mientras RF permanece inalterado; por ello el simulador también cuantifica cuándo conviene un respaldo híbrido VLC + RF.

---

## 2. Métricas Oficiales de Evaluación

En estricto cumplimiento del enunciado del proyecto (*"evaluar la eficiencia y la confiabilidad"*, límite de **hasta dos métricas**):

1. **Métrica 1 — Confiabilidad (`BER`):** Tasa de error de bit (post-FEC cuando se activa Reed-Solomon RS(255,239)). A partir del umbral industrial $\text{BER} \le 10^{-6}$ se derivan dos resúmenes operativos:
   - **Cobertura espacial [%]:** porcentaje del área libre de tránsito con $\text{BER} \le 10^{-6}$.
   - **Disponibilidad temporal [%]:** porcentaje de ubicaciones que mantienen $\text{BER} \le 10^{-6}$ en al menos el **99 %** de las realizaciones Monte Carlo (descartando colisiones sólidas físicamente imposibles donde un obstáculo móvil se superpone sobre la coordenada exacta del AGV).
2. **Métrica 2 — Eficiencia (`Goodput` y `Eficiencia Energética Incremental`):**
   - **Goodput efectivo libre de errores [Mb/s]:**
     $$G_{\text{eff}} = \frac{R_b}{10^6} \cdot R_c \cdot (1 - \text{PER}), \qquad \text{PER} = 1 - (1 - \text{BER})^{L_{\text{pkt}}}$$
     con $L_{\text{pkt}} = 1000\text{ bits}$ y tasa de código $R_c = 1$ (sin FEC) o $R_c = 239/255 \approx 0,9373$ (RS-FEC).
   - **Eficiencia Energética Incremental de Comunicación [Mbit/J]:**
     $$\eta_{\text{EE}} = \frac{G_{\text{eff}}}{P_{\text{comm,extra}}}$$
     donde $P_{\text{comm,extra}}$ contabiliza únicamente la potencia eléctrica adicional requerida para comunicar (drivers de modulación + receptor en VLC; AP industrial + radio RF en Wi-Fi).
3. **Variable física intermedia:** La **SNR de apertura de ojo [dB]** (peor caso de ISI, sección 3.4) se representa en mapas espaciales de diagnóstico.

> **Nota sobre la eficiencia.** $P_{\text{comm,extra}}$ es un **supuesto** (1 W por luminaria y 0,5 W de receptor en VLC; 12 W de AP y 1,5 W de radio en Wi-Fi). Las conclusiones de eficiencia deben leerse como condicionales a esos valores; el experimento E11 entrega el punto de equilibrio.

---

## 3. Ecuaciones de los Modelos Físicos Implementados

### 3.1 Iluminancia industrial (EN 12464-1) y canal LOS Lambertiano
Cada luminaria LED se modela como un emisor Lambertiano de orden $m = -\ln(2) / \ln(\cos \Phi_{1/2})$ (Komine & Nakagawa, 2004). La intensidad luminosa axial es $I_0 = \frac{m+1}{2\pi} P_{\text{opt}} \cdot \text{LER}$ y la iluminancia horizontal incluye tanto la componente directa (con bloqueo por máquinas cuando aplica) como las reflexiones difusas de primer orden en las paredes.

Para un fotodiodo de área $A_{\text{pd}}$, normal unitaria $\hat{n}_{\text{rx}}$, filtro $T_s$, concentrador no formador de imagen de ganancia $g(\psi) = n_c^2 / \sin^2(\Psi_c)$ (si $\psi \le \Psi_c$) y trayecto despejado de longitud $d$:

$$H_{\text{LOS}} = \frac{(m + 1) A_{\text{pd}}}{2\pi d^2} \cos^m(\phi) \, T_s \, g(\psi) \cos(\psi) \, \mathbb{I}(\psi \le \Psi_c) \, \mathbb{I}(\text{no bloqueado})$$

### 3.2 Reflexiones difusas de primer orden en paredes
Las cuatro paredes de la nave (reflectividad difusa $\rho = 0,4$, Uysal et al., 2017) se discretizan en parches de $dA = 0,5 \times 0,5\text{ m}$. Cada parche reirradia Lambertianamente hacia el receptor con ganancia diferencial:

$$dH^{(1)} = \frac{(m + 1) A_{\text{pd}}}{2\pi^2 d_1^2 d_2^2} \rho \, dA \, \cos^m(\phi) \cos(\alpha) \cos(\beta) \, T_s \, g(\psi) \cos(\psi)$$

llegando con retardo $(d_1 + d_2)/c$ y atenuación por polvo $e^{-\beta(d_1 + d_2)}$. En los escenarios S4 y S5 las reflexiones incluyen bloqueo tanto de maquinaria fija como de operarios y grúas móviles.

### 3.3 Polvo en suspensión: Extinción Beer-Lambert + Dispersión (*Scattering*) Henyey-Greenstein
A diferencia de un modelo de extinción pura (que produce la paradoja irreal de mejorar la SNR al aumentar el polvo a altas tasas), se combinan extinción balística y dispersión hacia adelante (*forward scattering*) fundamentada en Al-Kinani et al. (2018), Palacios Játiva et al. (2020, 2022) y Schneider et al. (2021):

1. **Coeficiente de extinción** (relación de Koschmieder a 550 nm): $\beta = 3,912 / V$, profundidad óptica $\tau_o = \beta d$.
2. **Componente balística:** $g_{\text{los}} = H_{\text{LOS}} e^{-\tau_o}$.
3. **Componente dispersada hacia adelante:** Con albedo de dispersión simple $\omega_0 = 0,65$ y factor de asimetría de Henyey-Greenstein $g_{\text{HG}} = 0,70$, la fracción angular que entra al cono del FOV es la integral analítica exacta:
   $$\Phi_{\text{HG}}(\Psi_c; g_{\text{HG}}) = \frac{1 - g_{\text{HG}}^2}{2 g_{\text{HG}}} \left(\frac{1}{1 - g_{\text{HG}}} - \frac{1}{\sqrt{1 + g_{\text{HG}}^2 - 2 g_{\text{HG}} \cos\Psi_c}}\right)$$
   $$g_{\text{scat}} = H_{\text{LOS}} \, \Phi_{\text{HG}}(\Psi_c) \, \left(1 - e^{-\omega_0 \tau_o}\right) e^{-(1 - \omega_0)\tau_o}$$
4. **Ensanchamiento temporal:** Los fotones dispersados llegan tras $t_{\text{LOS}}$ con retardo exponencial de escala $\tau_s = \frac{d}{c}(1 + \tau_o)\frac{1 - g_{\text{HG}}}{2}$, incrementando tanto el ruido *shot* como la ISI.

### 3.4 ISI por apertura de ojo con muestreo óptimo y ancho de banda del LED (revisión técnica 2)
Cada llegada $k$ (LOS de cada luminaria, cola de *scattering* representada por su instante medio $t_{\text{LOS}} + \tau_s$, y cada intervalo del histograma de reflexiones) tiene un tiempo $t_k$ y una ganancia $g_k$.

- **Instante de muestreo óptimo.** Un receptor integrar-y-descargar con ventana $[t_0, t_0 + T_b]$ recibe de una llegada desplazada $\Delta t_k = t_k - t_0$ la fracción $1 - |\Delta t_k|/T_b$ del bit actual y $|\Delta t_k|/T_b$ del bit vecino (anterior si $\Delta t_k > 0$, siguiente si $\Delta t_k < 0$). La amplitud total de ISI $\sum_k g_k |\Delta t_k|$ es mínima cuando $t_0$ es la **mediana ponderada** de los $t_k$.
- **Ancho de banda del LED blanco (Le Minh et al., 2009):** filtro paso bajo de 1.er orden con $\tau_{\text{LED}} = 1/(2\pi f_{3\text{dB}})$; la fracción útil es $\eta_{\text{LED}} = 1 - e^{-T_b/\tau_{\text{LED}}}$ y la cola $(1 - \eta_{\text{LED}})$ se suma a la ISI del bit anterior. $f_{3\text{dB}} = 3\text{ MHz}$ (LED fosforescente sin ecualizar), $15\text{ MHz}$ (pre-ecualización, diseño base) y $25\text{ MHz}$ (filtro azul, $T_{s,\text{azul}} = 0,55$ y $T_{s,\text{bg}} = 0,25$; valores supuestos).
- **BER con ISI determinista.** Con amplitudes de corriente $A_s$ (señal), $A_\ell$ (ISI tardía) y $A_e$ (ISI temprana) y ruido $\sigma^2 = \sigma^2_{\text{shot}} + \sigma^2_{\text{term}}$:
  $$\text{BER} = \frac{1}{4}\sum_{\pm,\pm} Q\!\left(\frac{A_s \pm A_\ell \pm A_e}{\sigma}\right), \qquad \text{SNR}_{\text{ojo}} = \frac{(A_s - A_\ell - A_e)^2}{\sigma^2}$$
- **Por qué cambió.** La versión anterior sincronizaba con la primera llegada y sumaba $I_{\text{ISI}}^2$ al ruido gaussiano (Komine & Nakagawa, 2004). Con 9 luminarias a distancias distintas la ISI resultaba ~1700 veces mayor que el ruido, la SNR quedaba fija (~26 dB) frente a temperatura y luz ambiente, y bloquear luminarias lejanas *mejoraba* la SNR. Ese modelo se conserva como `isi_model="gauss"` para comparación.

### 3.5 Suciedad depositada (LDD), temperatura de unión ($T_j$) y acoplamiento de $I_{\text{bg}}$
El fenómeno de *Luminaire Dirt Depreciation* (LDD) está documentado por Davis et al. (2019) en luminarias de exterior, que midieron depreciaciones más lentas que las usadas aquí. Las constantes siguientes y el acoplamiento térmico son **supuestos** para polvo industrial denso, evaluados con barridos en E5 y E6:
- **Transmitancia óptica:** $T_{\text{soil,rx}}(t) = e^{-t/\tau_{\text{rx}}}$ ($\tau_{\text{rx}} = 30\text{ d}$, lente orientado hacia arriba) y $T_{\text{soil,tx}}(t) = e^{-t/\tau_{\text{tx}}}$ ($\tau_{\text{tx}} = 120\text{ d}$, luminaria orientada hacia abajo).
- **Temperatura de unión del LED:**
  $$T_j = T_{\text{amb}} + \Delta T_{j,\text{nom}} \Big[1 + k_{\text{soil,th}}(1 - T_{\text{soil,tx}})\Big]$$
  con $\Delta T_{j,\text{nom}} = 30,25\,^\circ\text{C}$ ($0,55\,^\circ\text{C/W} \times 55\text{ W}$ térmicos), caída térmica del flujo $k_T = 0,4\,\%/^\circ\text{C}$ y protección térmica del driver sobre $T_{j,\text{derate}} = 85\,^\circ\text{C}$ ($1,5\,\%/^\circ\text{C}$).
- **Luz ambiente acoplada físicamente:**
  $$I_{\text{bg,eff}} = I_{\text{bg,0}} \cdot \left(\frac{A_{\text{rama}}}{1\text{ cm}^2}\right) \cdot T_{\text{soil,rx}}(t) \cdot T_{s,\text{filtro}} \cdot T_{s,\text{bg,azul}}$$

### 3.6 Ruido, SNR, Codificación RS(255,239) FEC y Receptor ADR Piramidal
- **Ruido *shot* + térmico** según Komine & Nakagawa (2004), incluyendo corriente oscura que se duplica cada $10\,^\circ\text{C}$.
- **Codificación FEC Reed-Solomon RS(255,239)** (IEEE Std 802.15.7-2018, $m=8\text{ bits/símbolo}$, $t=8\text{ símbolos}$, $R_c = 239/255$):
  $$\text{BER}_{\text{post}} = p_b \cdot I_{p_s}(t,\; n - t), \qquad p_s = 1 - (1 - p_b)^8$$
  donde $I_x(a, b)$ es la función beta incompleta regularizada.
- **Receptor con Diversidad Angular piramidal (ADR)** (Nuwanpriya et al., 2015): 5 fotodiodos (1 cenital + 4 inclinados a $35^\circ$ en azimutes $0^\circ, 90^\circ, 180^\circ, 270^\circ$) con FOV estrecho de $45^\circ$ y área por rama $A_{\text{pd}}/5$ (conservando el área total de silicio de $1\text{ cm}^2$), combinados mediante *Maximal Ratio Combining* con pesos $c_b = A_{s,b}/\sigma_b^2$; con el modelo de ojo la señal y la ISI combinadas son $\sum_b A_{s,b}^2/\sigma_b^2$ y $\sum_b A_{s,b}A_{\ell,b}/\sigma_b^2$ (y análogo para $A_e$), normalizadas por $\sqrt{\sum_b A_{s,b}^2/\sigma_b^2}$.

### 3.7 Referencia RF: Wi-Fi 2,4 GHz (3GPP TR 38.901 InF-SH) + EMI Middleton Class A
- **Pérdidas de propagación:** 3GPP TR 38.901 sec. 7.4.1 (*Indoor Factory – Sparse clutter, High BS*), determinando LOS/NLOS únicamente por maquinaria metálica pesada ($z \ge 2\text{ m}$) y modelando operarios interceptando el rayo con pérdida corporal de $4\text{ dB}$ (eliminando el doble conteo de obstrucciones).
- **Ruido impulsivo industrial Middleton Class A** (Middleton, 1977; Blackard et al., 1993):
  $$\text{BER}_{\text{ClassA}}(\gamma_0; A, \Gamma) = \sum_{m=0}^{M} e^{-A}\frac{A^m}{m!} \, Q\!\left(\sqrt{\frac{2\gamma_0}{1 + m/(A\Gamma)}}\right)$$
  donde $\gamma_0 = (E_b / N_{0,\text{th}})$ a la misma tasa objetivo $R_b = 10\text{ Mb/s}$ ($BW = 20\text{ MHz}$), $A = 0,1$ es el índice de impulsividad industrial y $\Gamma = 10^{-\text{EMI}_{\text{dB}}/10}$.
- **Mapa espacial de EMI:** Superposición en potencia de un piso de planta (`14 dB`) más la emisión local de horno de inducción, soldadura, prensa y CNCs con decaimiento $d^{-2}$.

---

## 4. Tabla de Trazabilidad Bibliográfica (Parámetro → Referencia)

| Parámetro / Modelo | Valor en simulador | Fuente científica verificada |
| --- | --- | --- |
| Canal LOS, 1.ª reflexión, $I_2, I_3, G_{\text{ol}}, \eta, \Gamma, g_m, R, I_{\text{bg}}$ | $0,562$; $0,0868$; $10$; $112\text{ pF/cm}^2$; $1,5$; $30\text{ mS}$; $0,53\text{ A/W}$; $1000/5100\text{ }\mu\text{A}$ | Komine & Nakagawa, *IEEE Trans. Consumer Electron.*, 2004 (`10.1109/TCE.2004.1277847`) |
| Reflectividad de paredes en celda industrial ($\rho$) | $0,4$ | Uysal, Miramirkhani et al., *IEEE Commun. Mag.*, 2017 (`10.1109/MCOM.2017.1600872CM`) |
| BW de LED blanco crudo / pre-ecualizado / con filtro azul | $3\text{ MHz}$ / $15\text{ MHz}$ / $25\text{ MHz}$ ($T_{s,\text{azul}} = 0,55$) | Le Minh et al., *IEEE Photon. Technol. Lett.*, 2009 (`10.1109/LPT.2009.2022413`) |
| Polvo (extinción Beer-Lambert + *scattering* HG $\omega_0, g_{\text{HG}}$) | $\omega_0 = 0,65$, $g_{\text{HG}} = 0,70$ | Al-Kinani et al., *IEEE COMST*, 2018 (`10.1109/COMST.2018.2838096`); Palacios Játiva et al., *IEEE Access* 2020 (`10.1109/ACCESS.2020.3030615`) & *MDPI Sensors* 2022 (`10.3390/s22072483`); Schneider et al., *Optics Express* 2021 (`10.1364/OE.421757`) |
| Depreciación por suciedad (LDD) y acople térmico $T_j$ | $\tau_{\text{rx}}=30\text{ d}$, $\tau_{\text{tx}}=120\text{ d}$, $k_{\text{soil,th}}=0,25$ (**supuestos**) | Concepto de LDD: Davis, Wilkerson & Kinzey, *LEUKOS*, 2019 (`10.1080/15502724.2018.1517598`), mediciones en exterior |
| Receptor ADR piramidal (5 ramas, inclinación y combinación MRC) | $\theta_{\text{tilt}}=35^\circ$, $\text{FOV}=45^\circ$ | Nuwanpriya, Ho & Chen, *IEEE JSAC*, 2015 (`10.1109/JSAC.2015.2432514`) |
| FEC Reed-Solomon RS(255,239) en VLC | $n=255, k=239, t=8$ | IEEE Std 802.15.7-2018 (`10.1109/IEEESTD.2019.8697198`) |
| Eficiencia energética verde y arquitectura VLC + Uplink IR | Uplink IR $850\text{ nm}$ ($200\text{ mW}$, $2\text{ Mb/s}$), Mbit/J | Jenila & Jeyachitra, *Digital Commun. & Networks* (ScienceDirect), 2021 (`10.1016/j.dcan.2020.09.004`); Matheus et al., *IEEE COMST*, 2019 (`10.1109/COMST.2019.2913348`) |
| Propagación RF fábrica (InF-SH) | $\sigma_{\text{LOS}}=4,3\text{ dB}$, $\sigma_{\text{NLOS}}=5,9\text{ dB}$ | 3GPP TR 38.901, sec. 7.4.1 |
| EMI impulsiva industrial en 2,4 GHz (Middleton Class A) | $A=0,1$; picos de horno/soldadura/VFD (**supuestos**) | Middleton, *IEEE Trans. EMC*, 1977 (`10.1109/TEMC.1977.303527`); Blackard, Rappaport & Bostian, *IEEE JSAC*, 1993 (`10.1109/49.233212`) |
| ISI por apertura de ojo y muestreo óptimo | mediana ponderada de llegadas; BER promediada sobre 4 patrones | Revisión técnica 2 (ver `docs/REVISION_TECNICA.md`) |
| Iluminación de áreas de trabajo industrial | $E_{\text{med}} \ge 300\text{ lx}$, $U_0 \ge 0,60$ | Norma EN 12464-1 |

---

## 5. Alcance y Limitaciones Declaradas

- Las reflexiones se modelan hasta primer orden sobre las cuatro paredes (en naves con superficies metálicas especulares pueden existir reflexiones de segundo orden adicionales en maquinaria).
- La dispersión por polvo se aproxima mediante un modelo semi-analítico de dispersión hacia adelante (*forward scattering* con función de fase de Henyey-Greenstein y cola exponencial) en lugar de trazado Monte Carlo fotón a fotón en 3D (como el usado en túneles mineros en Palacios Játiva et al., 2022), preservando tiempos de ejecución interactivos de pocos minutos.
- El traspaso (*handover*) vertical en la arquitectura híbrida se modela como selección ideal sin latencia de conmutación de capa MAC.
- En el modelo de ojo la ISI se agrupa en un bit anterior y uno siguiente (peor caso coherente); a tasas donde los retardos superan $T_b$ es una aproximación conservadora.
- No se modela el ruido de parpadeo (*flicker*) de luminarias ajenas que Schneider et al. (2021) miden en plantas reales; solo el ruido *shot* de la luz ambiente.
- En el simulador web (`docs/index.html`) los obstáculos móviles bloquean la línea de vista pero no las reflexiones.
- En el modo ADR se asume área total constante de silicio ($1\text{ cm}^2$ repartido en 5 ramas) y estimación ideal de relación señal a ruido por rama para MRC.
