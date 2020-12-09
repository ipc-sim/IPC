import pathlib
import argparse

import pymesh


def create_parser():
    parser = argparse.ArgumentParser(
        description="Convert a mesh to a IPC compatible tetrahedral mesh.")
    parser.add_argument(
        "input", metavar="path/to/input", type=pathlib.Path,
        help="path to input mesh")
    parser.add_argument(
        "--output", metavar="path/to/output.msh", type=pathlib.Path,
        default=None, help="path to output mesh")
    parser.add_argument(
        "--cellsize", type=float, default=1.0, dest="cell_size",
        help="max radius of the circumscribed sphere of the output tet")
    return parser


def parse_arguments():
    parser = create_parser()
    args = parser.parse_args()
    if args.output is None:
        args.output = args.input.with_suffix(".msh")
    return args


def main():
    args = parse_arguments()

    input_mesh = pymesh.load_mesh(str(args.input))
    if input_mesh.num_voxels == 0:  # tetrahedralize the mesh
        output_mesh = pymesh.tetrahedralize(input_mesh, args.cell_size)
    else:  # mesh is already a tet mesh, just convert to IPC compatible
        output_mesh = input_mesh

    # pymesh.save_mesh(str(args.output), output_mesh, ascii=True)
    with open(args.output, mode='w') as f:
        f.write("$MeshFormat\n4 0 8\n$EndMeshFormat\n")
        f.write("$Entities\n0 0 0 1\n")
        f.write("0 {:g} {:g} {:g} {:g} {:g} {:g} 0 0\n".format(
            *output_mesh.nodes.min(axis=0), *output_mesh.nodes.max(axis=0)))
        f.write("$EndEntities\n")
        f.write("$Nodes\n")
        f.write("1 {0:d}\n0 3 0 {0:d}\n".format(output_mesh.num_nodes))
        for i, node in enumerate(output_mesh.nodes):
            f.write("{:d} {:g} {:g} {:g}\n".format(i + 1, *node))
        f.write("$EndNodes\n")

        f.write("$Elements\n")
        f.write("1 {0:d}\n0 3 4 {0:d}\n".format(output_mesh.num_elements))
        f.write("\n")
        for i, element in enumerate(output_mesh.elements):
            f.write("{:d} {:d} {:d} {:d} {:d}\n".format(i + 1, *(element + 1)))
        f.write("$EndElements\n")

        f.write("$Surface\n")
        f.write("{:d}\n".format(output_mesh.num_faces))
        for face in output_mesh.faces:
            f.write("{:d} {:d} {:d}\n".format(*(face + 1)))
        f.write("$EndSurface\n")


if __name__ == "__main__":
    main()
