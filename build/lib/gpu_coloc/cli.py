import sys
from gpu_coloc import CLPA, coloc, format

def main():
    if "-r" in sys.argv or "--run" in sys.argv:
        sys.argv.remove("-r") if "-r" in sys.argv else sys.argv.remove("--run")
        coloc.main()
    if "-c" in sys.argv or "--clpa" in sys.argv:
        sys.argv.remove("-c") if "-c" in sys.argv else sys.argv.remove("--clpa")
        CLPA.main()
    elif "-f" in sys.argv or "--format" in sys.argv:
        sys.argv.remove("-f") if "-f" in sys.argv else sys.argv.remove("--format")
        format.main()
    else:
        print("Usage: gpu-coloc [-r|--run] or [-c|--clpa] or [-f|--format]")
        print("Use -r or --run to run the coloc script.")
        print("Use -c or --clpa to run the CLPA script.")
        print("Use -f or --format to run the format script.")
