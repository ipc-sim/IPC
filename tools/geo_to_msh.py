import sys
import json
import pathlib
import numpy
import pymesh

input_path = pathlib.Path(sys.argv[1])
output_path = input_path.with_suffix(".msh")

assert(input_path.exists())

with open(input_path, "r") as f:
    geo_mesh = json.load(f)

geo_dict = {}
assert(len(geo_mesh) % 2 == 0)
for i in range(0, len(geo_mesh), 2):
    geo_dict[geo_mesh[i]] = geo_mesh[i + 1]

attributes = {}
assert(len(geo_dict["attributes"]) % 2 == 0)
for i in range(0, len(geo_dict["attributes"]), 2):
    attributes[geo_dict["attributes"][i]] = geo_dict["attributes"][i + 1]

vertices = numpy.array(
    attributes["pointattributes"][0][1][7][5]).reshape(-1, 3)

faces = numpy.empty((0, 3))

voxels = numpy.array(geo_dict["topology"][1][1])
assert(voxels.size % 4 == 0)
voxels = voxels.reshape(-1, 4)

mesh = pymesh.form_mesh(vertices, faces, voxels)
pymesh.save_mesh(str(output_path), mesh, ascii=True)
