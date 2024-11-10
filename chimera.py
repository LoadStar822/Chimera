import argparse
import os
import shutil
import sys
import subprocess
from src.download import download
from src.profile import conversion2Krona
from src.profile import profile


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
    build_parser.add_argument("-m", "--mode", default="normal",
                              choices=["fast", "normal"],
                              help="Mode for building (choices: 'fast', 'normal'). Default is 'normal'.")
    build_parser.add_argument("-k", "--kmer", type=kmer_type, default=19,
                              help="Kmer size for building (must be between 1 and 31)")
    build_parser.add_argument("-w", "--window", type=int, default=31, help="Window size for building")
    build_parser.add_argument("-l", "--min-length", type=int, default=0, help="Minimum length sequence for building")
    build_parser.add_argument("-t", "--threads", type=int, default=32, help="Number of threads for building")
    build_parser.add_argument("--load-factor", type=float, default=0.95, help="Loading ratio of ICF")
    build_parser.add_argument("-M", "--max-hashes", type=int, default=2000000, help="Maximum number of hashes per taxid")
    build_parser.add_argument("-a", "--alpha", type=float, default=1.2, help="Alpha value for HICF")
    build_parser.add_argument("--relaxed-load-factor", type=float, default=0.95, help="Relaxed loading ratio of HICF")
    build_parser.add_argument("-f", "--filter", default="icf", choices=["icf", "hicf"],
                              help="Filter for building (choices: 'icf', 'hicf'). Default is 'icf'.")
    build_parser.add_argument("-c", "--fixed-cutoff", type=int, help="Fixed cutoff for minimizer (0 - 255)")
    build_parser.add_argument("-q", "--quiet", action="store_false", help="Quiet output")

    # Download and Build combined subcommand
    download_build_parser = subparsers.add_parser("download_and_build",
                                                  help="Download NCBI database sequences and resources and build a taxonomic sequence database")
    download_build_parser.add_argument("-o", "--output", default="ChimeraDB", help="Output file for building")
    download_build_parser.add_argument("-m", "--mode", default="normal",
                              choices=["fast", "normal"],
                              help="Mode for building (choices: 'fast', 'normal'). Default is 'normal'.")
    download_build_parser.add_argument("-k", "--kmer", type=kmer_type, default=19,
                                       help="Kmer size for building (must be between 1 and 31)")
    download_build_parser.add_argument("-w", "--window", type=int, default=31, help="Window size for building")
    download_build_parser.add_argument("-l", "--min-length", type=int, default=0,
                                       help="Minimum length sequence for building")
    download_build_parser.add_argument("-t", "--threads", type=int, default=32, help="Number of threads for building")
    download_build_parser.add_argument("--load-factor", type=float, default=0.95, help="Loading ratio of ICF")
    download_build_parser.add_argument("-M", "--max-hashes", type=int, default=2000000, help="Maximum number of hashes per taxid")
    download_build_parser.add_argument("-a", "--alpha", type=float, default=1.2, help="Alpha value for HICF")
    download_build_parser.add_argument("--relaxed-load-factor", type=float, default=0.95, help="Relaxed loading ratio of HICF")
    download_build_parser.add_argument("-f", "--filter", default="icf", choices=["icf", "hicf"],
                              help="Filter for building (choices: 'icf', 'hicf'). Default is 'icf'.")
    download_build_parser.add_argument("-c", "--fixed-cutoff", type=int,
                              help="Fixed cutoff for minimizer (0 - 255)")
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
        default="normal",
        choices=["fast", "normal"],
        help="Mode for classifying (choices: 'fast', 'normal'). Default is 'normal'."
    )
    classify_parser.add_argument("-f", "--filter", default="icf", choices=["icf", "hicf"],
                                 help="Filter for classifying (choices: 'icf', 'hicf'). Default is 'icf'.")
    classify_parser.add_argument("-b", "--batch-size", type=int, default=400, help="Batch size for classifying")
    group = classify_parser.add_mutually_exclusive_group()

    group.add_argument("-l", "--lca", action="store_true", help="Use LCA algorithm for classification")
    classify_parser.add_argument("-T", "--tax-file", help="Taxonomy file for LCA classification")

    group.add_argument("-e", "--em", action="store_true",  help="Use EM algorithm for classification")
    group.add_argument("-V", "--vem", action="store_true", default=True, help="Use VEM algorithm for classification (default)")
    classify_parser.add_argument("--em-iter", type=int, default=100, help="Number of EM iterations")
    classify_parser.add_argument("--em-threshold", type=float, default=0.001, help="EM threshold")

    group.add_argument("--none", action="store_true", help="Do not use LCA or EM for classification")

    classify_parser.add_argument("-q", "--quiet", action="store_false", help="Quiet output")

    # Profile subcommand
    profile_parser = subparsers.add_parser("profile", help="Generate sequence profile")
    profile_parser.add_argument("-i", "--input", nargs='+', required=True, help="Input file(s) for profiling")
    profile_parser.add_argument("-o", "--output", default="ChimeraProfile", help="Output file for profiling")
    profile_parser.add_argument("-k", "--krona", action="store_true", help="Generate Krona chart")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


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
        profile.process_file(args.input, args.output)
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
        command.extend(["-a", str(args.alpha)])
        command.extend(["--relaxed-load-factor", str(args.relaxed_load_factor)])
        command.extend(["-f", args.filter])
        if args.fixed_cutoff:
            command.extend(["-c", str(args.fixed_cutoff)])
        command.extend(["-M", str(args.max_hashes)])
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
        command.extend(["-f", args.filter])
        if args.lca:
            command.append("--lca")
            if not args.tax_file:
                raise ValueError("Taxonomy file must be provided when using LCA algorithm")
            command.extend(["--tax-file", args.tax_file])
        if args.em:
            command.append("-e")
            command.extend(["--em-iter", str(args.em_iter)])
            command.extend(["--em-threshold", str(args.em_threshold)])
        if args.vem:
            command.append("-V")
            command.extend(["--em-iter", str(args.em_iter)])
            command.extend(["--em-threshold", str(args.em_threshold)])
        if args.quiet:
            command.append("-q")

    # Execute the command using the provided run function
    download.run(command)


def main():
    chimera_path = get_chimera_path()
    args = parse_arguments()
    run_chimera(args, chimera_path)


if __name__ == "__main__":
    main()
