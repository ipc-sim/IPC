import numpy as np
from netCDF4 import Dataset


def get_matrix(var):
    shape = var.shape
    mat = np.zeros(shape)

    if len(shape) == 1:
        for i in range(shape[0]):
            mat[i] = var[i]
    elif len(shape) == 2:
        for i in range(shape[0]):
            for j in range(shape[1]):
                mat[i, j] = var[i, j]
    else:
        assert(False)

    return mat


def get_variable(data, name):
    if name in data.variables:
        return data.variables[name]

    return None


def get_variable_with_prefix(data, pre):
    res = []

    for v in data.variables:
        if v.startswith(pre):
            res.append(data.variables[v])

    return res


def read_coordinate(data, name):
    var = get_variable(data, name)

    if var is None:
        return None

    return get_matrix(var)


def read_points(data):
    var = get_variable(data, "coord")

    if var is not None:
        pts = get_matrix(var)

        return pts

    x = read_coordinate(data, "coordx")
    y = read_coordinate(data, "coordy")
    z = read_coordinate(data, "coordz")

    if x is None or y is None or z is None:
        return None

    pts = np.zeros((x.shape[0], 3))
    pts[:, 0] = x
    pts[:, 1] = y
    pts[:, 2] = z

    return pts


def read_connectivity(data):
    # static const std::string el = "elem_map",
    con = "connect"

    num_el_blk = data.dimensions["num_el_blk"].size
    if num_el_blk == 0:
        return None

    vars = get_variable_with_prefix(data, con)

    tetss = []

    for var in vars:
        elem_type_i = var.getncattr("elem_type")
        assert(elem_type_i == "TETRA")
        tets = get_matrix(var) - 1
        tetss.append(tets)

    s = 0

    for t in tetss:
        s += t.shape[0]

    res = np.zeros((s, 4), dtype=int)

    s = 0
    for t in tetss:
        res[s:s + t.shape[0], :] = t
        s += t.shape[0]

    return res


def writeMsh(nodes, elements, name):
    n_nodes = nodes.shape[0]
    n_els = elements.shape[0]
    assert(elements.shape[1] == 4)

    with open(name, 'w') as fid:
        fid.write('$MeshFormat\n')
        fid.write('2.2 0 8\n')
        fid.write('$EndMeshFormat\n')
        fid.write('$Nodes\n')
        fid.write('{}\n'.format(n_nodes))

        for i in range(n_nodes):
            fid.write('{} {} {} {}\n'.format(
                i + 1, nodes[i, 0], nodes[i, 1], nodes[i, 2]))

        fid.write('$EndNodes\n')
        fid.write('$Elements\n')
        fid.write('{}\n'.format(n_els))

        for i in range(n_els):
            fid.write('{} 4 2 0 0 {} {} {} {}\n'.format(
                i + 1, elements[i, 0] + 1, elements[i, 1] + 1, elements[i, 2] + 1, elements[i, 3] + 1))

        fid.write('$EndElements\n')


data = Dataset("utopiaComparison39K.e", "r+")

pts = read_points(data)
tets = read_connectivity(data)
writeMsh(pts, tets, "test.msh")

data.close()
