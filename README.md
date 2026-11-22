# Simulación de un sistema VLC en una nave industrial hostil

Proyecto 2 · Comunicaciones por luz visible para industrias inteligentes.

## Contenido

| Archivo | Qué hace |
| --- | --- |
| `vlc_core.py` | Modelos: canal Lambertiano LOS + reflexiones de primer orden, polvo (Beer-Lambert/Koschmieder), suciedad en ópticas, temperatura, obstáculos (segmento–caja), ruido shot/térmico (Komine 2004), ISI, BER OOK y referencia Wi-Fi 2,4 GHz (3GPP TR 38.901 InF-SH). |
| `run_all.py` | Ejecuta los 10 experimentos (E1–E10) y guarda figuras y CSV en `resultados/`. Tarda ≈ 2 min. Semilla fija: los resultados son reproducibles. |
| `demo_interactiva.py` | Demo en vivo para la presentación: mapa de SNR con controles de polvo, limpieza, temperatura, luz ambiente, operarios, tasa y FOV. |
| `resultados/` | Las 12 figuras, los CSV de cada barrido, `resumen_resultados.csv` y el registro de ejecución. |

## Cómo ejecutar

```bash
pip install -r requirements.txt
python run_all.py           # genera todos los resultados
python demo_interactiva.py  # abre la demo (requiere entorno gráfico)
```

## Experimentos

| Id | Experimento | Figura |
| --- | --- | --- |
| E1 | Planta de la nave e iluminancia | `E1_planta_iluminancia.png` |
| E2 | Caso de referencia: potencia, SNR, retardo RMS, respuesta impulsiva | `E2_*.png` |
| E3 | Seis escenarios (ideal, polvo, temperatura, luz ambiente, obstáculos, combinado) | `E3_*.png` |
| E4 | BER y cobertura vs visibilidad (10, 50, 100 Mb/s) | `E4_polvo_ber.png` |
| E5 | Cobertura vs días sin limpieza | `E5_suciedad_limpieza.png` |
| E6 | Efecto de la temperatura y componentes de ruido | `E6_temperatura.png` |
| E7 | Cobertura vs tasa de bits (ruido e ISI) | `E7_tasa_bits.png` |
| E8 | Monte Carlo de obstáculos móviles: FOV, nº de luminarias, Wi-Fi, híbrido | `E8_obstaculos_montecarlo.png` |
| E9 | VLC vs Wi-Fi frente a EMI y polvo | `E9_comparacion_wifi.png` |
| E10 | Diseño base vs mejorado en el escenario combinado | `E10_diseno_mejorado.png` |

## Métricas

1. **SNR** (dB).
2. **BER** de OOK-NRZ, BER = Q(√SNR). Umbral de confiabilidad: BER ≤ 10⁻⁶ (SNR ≥ 13,5 dB).

Cobertura = % del área con BER ≤ 10⁻⁶. Disponibilidad = % de realizaciones Monte Carlo con BER ≤ 10⁻⁶ en una ubicación (se considera confiable si es ≥ 99 %).

## Cambiar parámetros

Todos los parámetros están en la clase `Params` de `vlc_core.py` (unidades SI). Para un caso puntual:

```python
from vlc_core import *
p = replace(Params(), visibility=15.0, T_amb=60.0, Rb=20e6)
```

## Limitaciones del modelo

Solo extinción por polvo (sin luz dispersada), reflexiones de primer orden en paredes, sin ancho de banda del LED, EMI como ruido gaussiano equivalente, Wi-Fi y VLC sin codificación. Ver la sección 7 del informe.

## Referencias principales

- T. Komine, M. Nakagawa, IEEE Trans. Consumer Electronics, 50(1), 100–107, 2004.
- P. Palacios Játiva et al., IEEE Access, 8, 185445–185464, 2020.
- 3GPP TR 38.901, sección 7.4.1 (Indoor Factory).
