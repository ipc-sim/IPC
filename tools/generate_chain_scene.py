import sys

import numpy

scene_template = """\
energy NH
warmStart 0
size {:g}
time 10 5e-3
density 1e3
stiffness 1e6 0.4
script fallNoShift

shapes input {:d}
{:s}
selfCollisionOn

meshCO input/triMeshes/torus.obj  {:g} {:g} {:g}  1  50 0

constraintSolver interiorPoint

view orthographic
zoom 0.6

section SQP QP
energy FCR
time 10 1e-3
warmStart 1
constraintType graphics
constraintOffset 1e-3
tol 1
1e-3
section end
"""

if len(sys.argv) < 3:
    print("Usage: {} <num-links> <out-path>".format(sys.argv[0]))
    exit(1)

num_links = int(sys.argv[1])

delta_y = -0.51

shapes = "\n".join(
    ["input/tetMeshes/torus.msh  0 {:g} 0  0 {:g} 0  1 1 1".format(
        i * delta_y, 0 if i % 2 else 90) for i in range(num_links)])

scale = 1 + abs(delta_y) * (num_links - 1)

fixed_x = 0.5 if num_links > 1 else 0.190211 / 2
fixed_y = scale * 1 + 0.01
fixed_z = 0.5

scene = scene_template.format(
    scale, num_links, shapes, fixed_x, fixed_y, fixed_z)

with open(sys.argv[2], "w") as f:
    f.write(scene)
