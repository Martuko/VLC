# -*- coding: utf-8 -*-
"""
Construye la página del simulador 3D: inserta vlc_physics.js dentro de
index.template.html y escribe el resultado en ../docs/index.html
(la carpeta que GitHub Pages publica).

Uso (desde la raíz del repositorio):  python web/build.py
"""
import pathlib

here = pathlib.Path(__file__).resolve().parent
tpl = (here / "index.template.html").read_text(encoding="utf-8")
phys = (here / "vlc_physics.js").read_text(encoding="utf-8")
out = tpl.replace("/*__PHYSICS__*/", phys)
dest = here.parent / "docs" / "index.html"
dest.write_text(out, encoding="utf-8")
(here.parent / "docs" / ".nojekyll").write_text("", encoding="utf-8")
print(f"Página escrita en {dest} ({len(out.encode()) // 1024} KB)")
