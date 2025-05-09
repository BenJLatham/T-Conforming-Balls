import argparse
import gmsh
from T_Conf import Transfinite_Disk, Transfinite_Sphere


def main():
    parser = argparse.ArgumentParser(
        description="T-Conf: Symmetric transfinite mesh generator for disks and spheres using Gmsh"
    )

    subparsers = parser.add_subparsers(dest="command", help="Mesh geometry to generate")

    # Disk command
    parser_disk = subparsers.add_parser("disk", help="Generate a transfinite 2D disk mesh")
    parser_disk.add_argument("--a", type=float, required=True, help="Half-thickness of band")
    parser_disk.add_argument("--rInner", type=float, required=True, help="Radius of band center")
    parser_disk.add_argument("--rOuter", type=float, required=True, help="Outer radius")
    parser_disk.add_argument("--h_band", type=float, default=0.1, help="Band mesh size")
    parser_disk.add_argument("--h_outer", type=float, default=0.5, help="Outer mesh size")
    parser_disk.add_argument("--sectors", type=int, default=4, help="Number of angular sectors")
    parser_disk.add_argument("--out", type=str, default="disk.msh", help="Output mesh file")

    # Sphere command
    parser_sphere = subparsers.add_parser("sphere", help="Generate a transfinite 3D sphere mesh")
    parser_sphere.add_argument("--a", type=float, required=True, help="Half-thickness of band")
    parser_sphere.add_argument("--rInner", type=float, required=True, help="Radius of band center")
    parser_sphere.add_argument("--rOuter", type=float, required=True, help="Outer radius")
    parser_sphere.add_argument("--h_band", type=float, default=0.1, help="Band mesh size")
    parser_sphere.add_argument("--h_outer", type=float, default=0.1, help="Outer mesh size")
    parser_sphere.add_argument("--out", type=str, default="sphere.msh", help="Output mesh file")

    args = parser.parse_args()

    if args.command == "disk":
        Transfinite_Disk(
            a=args.a,
            rInner=args.rInner,
            rOuter=args.rOuter,
            h_band=args.h_band,
            h_outer=args.h_outer,
            n_sectors=args.sectors,
        )
        gmsh.model.mesh.generate(2)
        gmsh.write(args.out)
        gmsh.finalize()
        print(f"Disk mesh written to {args.out}")

    elif args.command == "sphere":
        Transfinite_Sphere(
            a=args.a,
            rInner=args.rInner,
            rOuter=args.rOuter,
            h_band=args.h_band,
            h_outer=args.h_outer,
        )
        gmsh.model.mesh.generate(3)
        gmsh.write(args.out)
        gmsh.finalize()
        print(f"Sphere mesh written to {args.out}")

    else:
        parser.print_help()


if __name__ == "__main__":
    main()
