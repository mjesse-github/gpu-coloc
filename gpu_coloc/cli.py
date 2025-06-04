import sys
from gpu_coloc import coloc, format

def main():
    if "-r" in sys.argv or "--run" in sys.argv:
        sys.argv.remove("-r") if "-r" in sys.argv else sys.argv.remove("--run")
        coloc.main()
    elif "-f" in sys.argv or "--format" in sys.argv:
        sys.argv.remove("-f") if "-f" in sys.argv else sys.argv.remove("--format")
        format.main()
    else:
        print("Usage: gpu-coloc [-r|--run] or [-f|--format]")
        print("Use -r or --run to run the coloc script.")
        print("Use -f or --format to run the format script.")
