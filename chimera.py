import argparse
import os
import shutil
import subprocess
import src.download.download as download
import src.profile.conversion2Krona as conversion2Krona


def get_chimera_path():
    # Use which command to find the chimera program path
    chimera_path = shutil.which("Chimera")
    if chimera_path is None:
        raise FileNotFoundError("Cannot find 'Chimera' executable. Please ensure it is installed and in your PATH.")
    return chimera_path


def kmer_type(x):
    x = int(x)
    if x < 1 or x > 31:
        raise argparse.ArgumentTypeError("Kmer size must be between 1 and 31")
    return x


def parse_arguments():
    parser = argparse.ArgumentParser(description="Chimera - A versatile tool for metagenomic classification")

    parser.add_argument("-v", "--version", action="store_true", help="Show version information")

    subparsers = parser.add_subparsers(dest="command", required=False)

    # Download subcommand
    download_parser = subparsers.add_parser("download", help="Download NCBI database sequences and resources")

    # Build subcommand
    build_parser = subparsers.add_parser("build", help="Build a taxonomic sequence database")
    build_parser.add_argument("-i", "--input", required=True, help="Input file for building")
    build_parser.add_argument("-o", "--output", default="ChimeraDB", help="Output file for building")
    build_parser.add_argument("-m", "--mode", default="default", help="Mode for building")
    build_parser.add_argument("-k", "--kmer", type=kmer_type, default=19,
                              help="Kmer size for building (must be between 1 and 31)")
    build_parser.add_argument("-w", "--window", type=int, default=31, help="Window size for building")
    build_parser.add_argument("-l", "--min-length", type=int, default=0, help="Minimum length sequence for building")
    build_parser.add_argument("-t", "--threads", type=int, default=32, help="Number of threads for building")
    build_parser.add_argument("--load-factor", type=float, default=0.95, help="Loading ratio of ICF")
    build_parser.add_argument("-q", "--quiet", action="store_false", help="Quiet output")

    # Download and Build combined subcommand
    download_build_parser = subparsers.add_parser("download_and_build",
                                                  help="Download NCBI database sequences and resources and build a taxonomic sequence database")
    download_build_parser.add_argument("-o", "--output", default="ChimeraDB", help="Output file for building")
    download_build_parser.add_argument("-m", "--mode", default="default", help="Mode for building")
    download_build_parser.add_argument("-k", "--kmer", type=kmer_type, default=19,
                                       help="Kmer size for building (must be between 1 and 31)")
    download_build_parser.add_argument("-w", "--window", type=int, default=31, help="Window size for building")
    download_build_parser.add_argument("-l", "--min-length", type=int, default=0,
                                       help="Minimum length sequence for building")
    download_build_parser.add_argument("-t", "--threads", type=int, default=32, help="Number of threads for building")
    download_build_parser.add_argument("--load-factor", type=float, default=0.95, help="Loading ratio of ICF")
    download_build_parser.add_argument("-q", "--quiet", action="store_false", help="Quiet output")

    # Classify subcommand
    classify_parser = subparsers.add_parser("classify", help="Classify sequences")
    classify_parser.add_argument("-i", "--single", nargs='+',
                                 help="Input files for classifying (supports multiple files)")
    classify_parser.add_argument("-p", "--paired", nargs='+', help="Paired input files for classifying")
    classify_parser.add_argument("-o", "--output", default="ChimeraClassify", help="Output file for classifying")
    classify_parser.add_argument("-d", "--database", required=True, help="Database file for classifying")
    classify_parser.add_argument("-s", "--shot-threshold", type=float, default=0.7,
                                 help="Shot threshold for classifying")
    classify_parser.add_argument("-t", "--threads", type=int, default=32, help="Number of threads for classifying")
    classify_parser.add_argument(
        "-m",
        "--mode",
        default="fast",
        choices=["fast", "normal"],
        help="Mode for classifying (choices: 'fast', 'normal'). Default is 'fast'."
    )
    classify_parser.add_argument("-b", "--batch-size", type=int, default=400, help="Batch size for classifying")
    classify_parser.add_argument("-q", "--quiet", action="store_false", help="Quiet output")

    # Profile subcommand
    profile_parser = subparsers.add_parser("profile", help="Generate sequence profile")
    profile_parser.add_argument("-i", "--input", nargs='+', required=True, help="Input file for profiling")
    profile_parser.add_argument("-o", "--output", default="ChimeraProfile", help="Output file for profiling")
    profile_parser.add_argument("-k", "--krona", action="store_true", help="Generate Krona chart")

    return parser.parse_args()


def validate_arguments(args):
    if not args.version and args.command is None:
        raise ValueError("You must specify one of the following commands: download, build, classify, profile")


def run_chimera(args, chimera_path):
    if args.command == "download":
        download.download(interactive=True)
        return 0

    if args.command == "download_and_build":
        options = download.download(interactive=True)
        args.command = "build"
        args.input = os.path.join(options.output_dir, "target.tsv")

    if args.command == "profile":
        if args.krona:
            print("Converting file to Krona format...")
            conversion2Krona.convert_multiple_files_to_krona_format(args.input, args.output)
            print("Conversion completed.")
            print("Generating Krona chart...")
            download.run(["ktImportText", args.output + ".tsv", "-o", args.output + ".html"])
            print("Krona chart generated.")
        return 0


    command = [chimera_path]

    if args.version:
        command.append("--version")

    # Add subcommand
    if args.command:
        command.append(args.command)

    # Add options and parameters based on parsed arguments
    if args.command == "download":
        if args.build:
            command.append("-b")

    if args.command == "build":
        command.extend(["-i", args.input])
        command.extend(["-o", args.output])
        command.extend(["-m", args.mode])
        command.extend(["-k", str(args.kmer)])
        command.extend(["-w", str(args.window)])
        command.extend(["-l", str(args.min_length)])
        command.extend(["-t", str(args.threads)])
        command.extend(["--load-factor", str(args.load_factor)])
        if args.quiet:
            command.append("-q")

    elif args.command == "classify":

        if args.single:
            command.extend(["-i"] + args.single)
        if args.paired:
            command.extend(["-p"] + args.paired)
        command.extend(["-o", args.output])
        command.extend(["-d", args.database])
        command.extend(["-s", str(args.shot_threshold)])
        command.extend(["-t", str(args.threads)])
        command.extend(["-m", args.mode])
        command.extend(["-b", str(args.batch_size)])
        if args.quiet:
            command.append("-q")

    # Execute the command using the provided run function
    download.run(command)


if __name__ == "__main__":
    try:
        chimera_path = get_chimera_path()
        args = parse_arguments()
        validate_arguments(args)
        run_chimera(args, chimera_path)
    except ValueError as e:
        print(f"Error: {e}")
    except FileNotFoundError as e:
        print(f"Error: {e}")
