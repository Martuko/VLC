# Revisión técnica 2: hallazgos, correcciones y simulador 3D

Esta revisión parte de la versión mejorada del repositorio (la de `REGISTRO_DE_MEJORAS.md`) y responde a tres preguntas: si los resultados son reproducibles, si los modelos y las referencias sostienen las conclusiones, y cómo llevar el simulador a una página web en 3D.

**Resumen.** La versión mejorada es reproducible y aporta modelos valiosos, pero su modelo de interferencia intersimbólica (ISI) dejaba el enlace saturado por ISI. Por eso la temperatura y la luz ambiente casi no tenían efecto, y bloquear trayectos con obstáculos *mejoraba* la SNR. Se corrigió con un modelo de apertura de ojo, se reemplazó una referencia no verificable, se precisó otra mal aplicada y se agregó un análisis de sensibilidad de la eficiencia (E11). Con el modelo corregido VLC resulta más robusto de lo reportado antes y el receptor ADR deja de ser imprescindible: un fotodiodo único con FOV de 80° logra la misma disponibilidad.

---

## 1. Reproducibilidad

| Verificación | Resultado |
| --- | --- |
| `pytest tests/` de la versión mejorada | 13 de 13 aprobadas |
| `python run_all.py` frente a `log_ejecucion.txt` entregado | Idéntico número a número (solo difieren los saltos de línea CRLF/LF) |
| Tiempo de ejecución | ≈ 2 min (ahora ≈ 3 min por el modelo de ojo y E11) |

## 2. Hallazgo principal: enlace saturado por ISI

### 2.1 Evidencia

En el modelo anterior la ISI se trataba como ruido gaussiano, $\text{SNR} = I_{\text{sig}}^2 / (\sigma^2_{\text{shot}} + \sigma^2_{\text{term}} + I_{\text{ISI}}^2)$, y el receptor se sincronizaba con la **primera** llegada. Con 9 luminarias a distancias distintas, la luz de las luminarias lejanas llega 10–35 ns después y se contaba como ISI. Desglose en el centro de la nave a 10 Mb/s:

| Caso | Ruido shot [A²] | Ruido térmico [A²] | ISI² [A²] | SNR (modelo anterior) |
| --- | --- | --- | --- | --- |
| 25 °C | 2,5·10⁻¹⁵ | 1,0·10⁻¹⁶ | 4,4·10⁻¹² | 26,1 dB |
| 70 °C | 2,3·10⁻¹⁵ | 1,2·10⁻¹⁶ | 1,8·10⁻¹² | 26,1 dB |
| Sol directo (5100 µA) | 9,9·10⁻¹⁵ | 1,0·10⁻¹⁶ | 4,4·10⁻¹² | 26,1 dB |

La ISI era unas 1700 veces mayor que todo el ruido, así que la SNR quedaba fija en la razón señal/ISI. Consecuencias:

1. **E6 (temperatura)** reportaba −0,33 dB a 70 °C aunque el LED perdía 37 % de potencia: era un artefacto.
2. **S3 (luz ambiente)** no tenía efecto, en contradicción con Schneider et al. (2021), que el propio repositorio cita y que identifica la interferencia de fuentes de luz como uno de los efectos más fuertes en planta.
3. **Paradoja de obstáculos:** en el escenario S4, la SNR con obstáculos era *mayor* que sin ellos en el 65 % de los puntos libres (mediana 30,1 dB frente a 28,7 dB en la misma área), porque bloquear luminarias lejanas eliminaba "ISI". Es la misma clase de paradoja que el registro anterior dice haber eliminado para el polvo.

### 2.2 Corrección: modelo de apertura de ojo (`Params.isi_model = "eye"`)

1. **Instante de muestreo óptimo.** Un receptor integrar-y-descargar con ventana $[t_0, t_0+T_b]$ recibe de una llegada desplazada $\Delta t$ la fracción $1-|\Delta t|/T_b$ del bit actual y $|\Delta t|/T_b$ del vecino (anterior si llega tarde, siguiente si llega antes). La amplitud total de ISI $\sum g_k |\Delta t_k|$ se minimiza con $t_0$ igual a la **mediana ponderada** de los tiempos de llegada.
2. **ISI determinista, no gaussiana.** Con amplitudes de señal $A_s$, ISI tardía $A_\ell$ (incluye la cola del LED) e ISI temprana $A_e$, la BER de OOK es el promedio sobre los cuatro patrones de bits vecinos:

$$\text{BER} = \frac{1}{4}\sum_{\pm,\pm} Q\!\left(\frac{A_s \pm A_\ell \pm A_e}{\sigma}\right)$$

3. **SNR reportada:** la de peor caso (apertura de ojo), $\text{SNR}_{\text{ojo}} = (A_s - A_\ell - A_e)^2/\sigma^2$; también se entrega `snr_noise_db` (sin ISI) para diagnóstico.
4. **ADR con MRC:** pesos $c_b = A_{s,b}/\sigma_b^2$ por rama; la ISI combinada se calcula con los mismos pesos.
5. El modelo anterior se conserva con `isi_model="gauss"` para comparar.

Efecto en el centro de la nave: 58,1 dB a 25 °C, 54,5 dB a 70 °C y 52,3 dB con sol directo. Puntos donde los obstáculos mejoran la SNR: 0 %.

### 2.3 Tests nuevos (16 en total)

- `test_eye_ber_limits`: sin ISI la BER es $Q(A/\sigma)$; con el ojo cerrado no baja de 1/4.
- `test_obstacles_never_improve_snr`: regresión de la paradoja de obstáculos.
- `test_snr_responds_to_noise_sources`: el calor y la luz ambiente reducen la SNR (y el modelo gaussiano conservado reproduce la saturación).

## 3. Verificación de referencias

| Referencia | Estado | Acción |
| --- | --- | --- |
| Komine & Nakagawa, IEEE TCE 2004 | Verificada | Ninguna |
| Schneider et al., Optics Express 29(11), 2021 | Verificada | Ninguna; ahora el modelo es coherente con su hallazgo sobre luz ambiente |
| Palacios Játiva et al., IEEE Access 2020 y Sensors 2022 | Verificadas | Ninguna |
| Jenila & Jeyachitra, Digital Commun. Networks 7(3), 2021 | Verificada | Ninguna |
| Davis, Wilkerson & Kinzey, LEUKOS 15(1), 2019 | Existe, pero mide suciedad en luminarias **de exterior**; no trata la temperatura de unión | Se cita solo para el concepto de LDD; $\tau_{rx}$, $\tau_{tx}$ y $k_{\text{soil,th}}$ se declaran supuestos |
| Oyedum et al., IEEE GLOBECOM 2025 (DOI `…11432514`) | **No verificable**: solo se encontró un proyecto de verano 2026 del WINLAB con ese autor | Eliminada; reemplazada por Blackard, Rappaport & Bostian, IEEE JSAC 11(7), 1993 (`10.1109/49.233212`) |
| Valores pico de EMI de `EMI_SOURCES` | Atribuidos a "mediciones NIST/Oyedum" sin respaldo | Declarados como supuestos de diseño |
| Uysal 2017, Le Minh 2009, Nuwanpriya 2015, Al-Kinani 2018, Matheus 2019, Middleton 1977, 3GPP TR 38.901, IEEE 802.15.7 | Publicaciones conocidas y coherentes con su uso | Ninguna |

## 4. La eficiencia (Mbit/J) depende de supuestos: experimento E11

La métrica divide el goodput por una potencia incremental **supuesta**: 1 W por luminaria para el driver VLC (repartido como 0,56 W en la grilla 4×4) frente a 12 W del punto de acceso Wi-Fi. E11 calcula el consumo por luminaria a partir del cual VLC deja de ser más eficiente:

| AP Wi-Fi | Eficiencia Wi-Fi + FEC | Equilibrio, 9 luminarias | Equilibrio, 16 luminarias ADR |
| --- | --- | --- | --- |
| 6 W | 1,24 Mbit/J | 0,84 W/luminaria | 0,42 W/luminaria |
| 12 W | 0,69 Mbit/J | 1,56 W/luminaria | 0,80 W/luminaria |
| 15 W | 0,56 Mbit/J | 1,91 W/luminaria | 0,99 W/luminaria |

La ventaja de eficiencia de VLC es real solo si el driver de modulación consume menos que esos valores. Así debe presentarse: como condición, no como resultado de la simulación.

## 5. Afirmaciones del README que se precisaron

- "La grilla 4×4 cumple EN 12464-1": cierto en la nave vacía (302 lx, $U_0$ = 0,61). Con las máquinas, la uniformidad en la zona de tránsito es 0,05, peor que la grilla 3×3 (0,08).
- "VLC 100 % con ADR + 4×4": es el resultado de E9 (aire limpio, 25 °C). En el escenario combinado hostil de E10 la disponibilidad es 78,9 %.
- "El ADR es la solución": con el modelo corregido, un fotodiodo único con FOV de 80° alcanza el 100 % en E8, y en E10 el ADR+4×4 (78,9 %) no supera al fotodiodo único con FEC (79,4 %).

## 6. Resultados antes y después de la corrección

| Indicador | Versión mejorada | Revisión 2 |
| --- | --- | --- |
| E2 · SNR mediana a 10 Mb/s | 28,7 dB | 56,2 dB |
| E3 · S2 calor 70 °C (SNR mediana) | 28,6 dB (sin efecto) | 52,5 dB (−3,7 dB) |
| E3 · S3 sol directo (SNR mediana) | 28,6 dB (sin efecto) | 50,2 dB (−6,0 dB) |
| E3 · S4 obstáculos (cobertura) | 92,6 % | 95,8 % |
| E3 · S5 combinado (cobertura) | 87,0 % | 89,7 % |
| E4 · polvo V = 5 m a 10 Mb/s (cobertura) | 59,1 % | 100 % |
| E6 · ΔSNR a 70 °C / 85 °C | −0,33 / −1,01 dB | −4,26 / −8,03 dB |
| E8 · 40 operarios, 1 fotodiodo FOV 60° | 70,0 % | 95,6 % |
| E8 · 40 operarios, 1 fotodiodo FOV 80° | 62,8 % | 100 % |
| E9 · VLC base 9 luminarias (EMI espacial) | 84,4 % | 99,4 % |
| E10 · Base / +Filtro / +FEC / ADR+4×4 | 60,0 / 62,2 / 71,7 / 78,9 % | 62,8 / 69,4 / 79,4 / 78,9 % |
| E10 · Wi-Fi + FEC / Híbrido | 70,6 / 99,4 % | 70,6 / 99,4 % (sin cambio: no depende del modelo VLC) |

**Nueva lectura de los resultados.** El factor dominante para VLC es el **bloqueo de línea de vista** (sombras detrás de estanterías y máquinas) y, en segundo lugar, la **suciedad acumulada** (con V = 10 m la cobertura cae al 55 % a los 120 días). El polvo en suspensión, el calor y la luz solar consumen margen (5–15 dB) pero no cortan el enlace a 10 Mb/s. Wi-Fi sigue limitado por la EMI cerca de horno, soldadura y variadores. La arquitectura híbrida mantiene el 99,4 %.

## 7. Simulador 3D web (`docs/index.html`)

- **Física:** `web/vlc_physics.js` es un port directo de `vlc_core.py` con el modelo de ojo (canal LOS con bloqueo, reflexiones de primer orden, polvo con scattering, suciedad, temperatura de unión, ruido, ancho de banda del LED, RS-FEC, ADR con MRC, Wi-Fi 3GPP InF-SH + Middleton Class A).
- **Validación:** `node web/validar_contra_python.js` compara con `web/referencia_python.json` (generado por `python web/generar_referencia.py`). En seis escenarios las diferencias son menores a 0,05 dB en la mediana y nulas en cobertura; Wi-Fi coincide dentro del ruido Monte Carlo.
- **Simplificación declarada:** en la web, operarios y grúas bloquean la línea de vista pero no las reflexiones en paredes (en Python sí, en S4/S5).
- **Construcción:** `python web/build.py` inserta la física en `web/index.template.html` y escribe `docs/index.html` (autocontenido; Three.js r128 y fuentes desde CDN).
- **Despliegue:** GitHub Pages con origen `main` → carpeta `/docs` (ver README).

## 8. Recomendaciones pendientes

1. **Ruido de parpadeo de luminarias ajenas.** Schneider et al. (2021) miden ruido de *flicker* de fuentes artificiales en planta; el modelo solo incluye ruido shot DC de la luz ambiente.
2. **Reflexiones en máquinas metálicas** (hoy solo paredes, primer orden).
3. **Ablación de un factor a la vez en E10:** el paso 2 cambia filtro y limpieza simultáneamente.
4. **Escalado del consumo en la grilla 4×4:** repartir 9 W entre 16 drivers favorece a la 4×4; lo más conservador es 1 W por luminaria (ver E11).
5. **`PyQt6` en `requirements.txt`:** solo hace falta si el backend Tk no está disponible; pesa más de 100 MB.

## 9. Preguntas probables en la defensa

- *¿Por qué cambió tanto la SNR entre versiones?* Porque la versión anterior trataba la luz de las luminarias lejanas como ruido aleatorio y sincronizaba con la primera llegada. El modelo de ojo muestrea en el instante óptimo y trata la ISI como lo que es: un efecto determinista de los bits vecinos.
- *¿Entonces el polvo ya no importa?* Importa: resta hasta 15 dB a V = 5 m. Pero a 10 Mb/s el margen es de más de 40 dB, así que no corta el enlace; a tasas más altas o con ópticas sucias sí.
- *¿Por qué no usar ADR?* Con este modelo, un fotodiodo con FOV de 80° ve más luminarias y logra la misma disponibilidad con un receptor más simple. El ADR conviene si se necesita FOV estrecho (por ejemplo, para rechazar luz solar lateral).
- *¿VLC es más eficiente que Wi-Fi?* Depende del consumo del driver: ver la tabla de E11.
