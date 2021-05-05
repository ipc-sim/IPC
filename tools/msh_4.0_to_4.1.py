import pathlib
import argparse
from collections import OrderedDict


def read_msh40(msh_path):
    V = OrderedDict()  # vertices
    T = OrderedDict()  # tets
    F = []

    with open(msh_path, mode='r') as f:
        lines = f.read().split("\n")

    i = 0
    while i < len(lines):
        if lines[i] == "$MeshFormat":
            version, is_binary, data_size = lines[i + 1].split()
            if version != "4" and version != "4.0":
                raise Exception(f"Invalid input MSH version: {version}")
            if int(is_binary):
                raise Exception(
                    f"Only able to convert ASCII format MSH files!")
            i += 3
        elif lines[i] == "$Entities":
            print("Skipping $Entities section")
            i += 4
        elif lines[i] == "$Nodes":
            _, num_nodes = lines[i + 1].split()
            _, _, _, num_nodes = lines[i + 2].split()
            num_nodes = int(num_nodes)
            i += 3
            for j in range(num_nodes):
                id, *node = lines[i + j].split()
                V[id] = node
            i += num_nodes + 1
            pass
        elif lines[i] == "$Elements":
            _, num_tets = lines[i + 1].split()
            _, _, _, num_tets = lines[i + 2].split()
            num_tets = int(num_tets)
            i += 3
            for j in range(num_tets):
                id, *tet = lines[i + j].split()
                T[id] = tet
            i += num_tets + 1
        elif lines[i] == "$Surface":
            num_tris = int(lines[i + 1])
            for j in range(num_tris):
                F.append(lines[i + j + 2].split())
            i += num_tris + 3
        elif lines[i] != "":
            print(f"Unkown line: {lines[i]}")
            i += 1
        else:
            i += 1

    return V, T, F


def write_msh41(outpath, V, T, F=None):
    with open(outpath, mode='w') as f:
        f.write("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")

        f.write("$Nodes\n")
        f.write("1 {:d} {} {}\n".format(
            len(V),
            min(V.keys(), key=lambda x: int(x)),
            max(V.keys(), key=lambda x: int(x))))
        f.write("3 0 0 {:d}\n".format(len(V)))
        f.write("{}\n".format('\n'.join(V.keys())))
        for node in V.values():
            f.write("{} {} {}\n".format(*node))
        f.write("$EndNodes\n")

        f.write("$Elements\n")
        f.write("1 {:d} {} {}\n".format(
            len(T),
            min(T.keys(), key=lambda x: int(x)),
            max(T.keys(), key=lambda x: int(x))))
        f.write("3 0 4 {}\n".format(len(T)))
        for id, elm in T.items():
            f.write("{} {} {} {} {}\n".format(id, *elm))
        f.write("$EndElements\n")

        if F:
            f.write("$Surface\n")
            f.write(f"{len(F):d}\n")
            for face in F:
                f.write("{} {} {}\n".format(*face))
            f.write("$EndSurface\n")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Convert a MSH 4.0 to the 4.1 specification.")
    parser.add_argument(
        "-i", "--input", metavar="path/to/input", type=pathlib.Path,
        dest="input", help="path to input MSH 4.0 file(s)", nargs="+",
        required=True)
    parser.add_argument(
        "-o", "--output", metavar="path/to/output/", type=pathlib.Path,
        dest="output", default=pathlib.Path("output"),
        help="path to output directory")
    parser.add_argument(
        "--keep-surface", action="store_true", default=False,
        help="keep the $Surface section")
    return parser.parse_args()


def main():
    args = parse_args()

    tet_mesh_dir = (
        pathlib.Path(__file__).resolve().parents[1] / "input" / "tetMeshes")

    for input in args.input:
        if args.output.is_dir():
            try:
                output = input.resolve().relative_to(tet_mesh_dir)
            except:
                output = input.name
            output = args.output / output
        else:
            output = args.output
        output.parent.mkdir(parents=True, exist_ok=True)

        print(f"Converting {input} to {output}")
        try:
            write_msh41(output,
                        *(read_msh40(input)[:3 if args.keep_surface else 2]))
        except Exception as err:
            print(err)


if __name__ == "__main__":
    main()
