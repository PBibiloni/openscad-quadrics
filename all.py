import os
import subprocess
from pathlib import Path

import cone, ellipsoid, hyperbolicparaboloid, hyperboloid, paraboliccylinder

OPENSCAD_EXEC = r'C:\Program Files\OpenSCAD\openscad.exe'

for w in [1.6, 2, 2.8]:
    cone.generate(0.5, 0.5, width_mm=w)
    cone.generate(0.5, 1, width_mm=w)
    cone.generate(2, 1, width_mm=w)

    ellipsoid.generate(width_mm=w)
    hyperbolicparaboloid.generate(width_mm=w)
    hyperboloid.generate(width_mm=w)
    paraboliccylinder.generate(width_mm=w)


root = Path(__file__).parent
for path_scad_dir in [p for p in root.iterdir() if p.name.startswith('scad_')]:
    path_stl_dir = root / path_scad_dir.name.replace('scad_', 'stl_')
    os.makedirs(path_stl_dir, exist_ok=True)
    for path_scad_script in [p for p in path_scad_dir.iterdir() if p.suffix == '.scad']:
        for out_ext in ['.stl', '.png']:
            path_out = path_stl_dir / path_scad_script.name.replace('.scad', out_ext)
            cmd = f'"{OPENSCAD_EXEC}" "{path_scad_script}" --imgsize=1024,1024 --render -o "{path_out}"'
            print(f'EXECUTING: {cmd}')
            subprocess.run(cmd, shell=True)
