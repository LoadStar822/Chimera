import argparse
import re
import subprocess
import sys
import time
import os
import shlex
from pathlib import Path
import pandas as pd
from multitax import NcbiTx

# Constants definition
VALID_DATABASES = ["genbank", "refseq"]
VALID_ORGANISM_GROUPS = [
    "archaea", "bacteria", "fungi", "human", "invertebrate", "metagenomes",
    "other", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral"
]
VALID_ASSEMBLY_LEVELS = ["complete genome", "chromosome", "scaffold", "contig"]
VALID_REFSEQ_CATEGORIES = ["reference genome", "na"]
VALID_FILE_TYPES = ["genomic.fna.gz", "assembly_report.txt", "protein.faa.gz", "genomic.gbff.gz"]
VALID_TAXONOMY_MODES = ["ncbi", "gtdb"]
VALID_DOWNLOADERS = ["wget", "curl"]

DEFAULT_DATABASE = "refseq"
DEFAULT_ASSEMBLY_LEVEL = "complete genome"
DEFAULT_REFSEQ_CATEGORY = "reference genome"
DEFAULT_FILE_TYPE = "genomic.fna.gz"
DEFAULT_OUTPUT_DIR = "./genome_output"
DEFAULT_THREADS = "1"
DEFAULT_TAXONOMY_MODE = "ncbi"
DEFAULT_RETRY_ATTEMPTS = "3"
DEFAULT_DOWNLOADER = "wget"


def validate_input(prompt, valid_options, default=None, allow_empty=False):
    """
    Validate user input
    
    Parameters:
    - prompt (str): Message to prompt user
    - valid_options (list): List of valid options
    - default (str): Default value
    - allow_empty (bool): Whether to allow empty input
    
    Returns:
    - str: Validated input value
    """
    while True:
        user_input = input(prompt).strip()
        if not user_input and default is not None:
            return default
        if allow_empty and not user_input:
            return ""
        entries = [entry.strip() for entry in user_input.split(",")]
        if all(entry in valid_options for entry in entries):
            return ",".join(entries)
        print(f"\nInvalid input. Please choose from: {', '.join(valid_options)}")


def prompt_user(options):
    """
    Interactively prompt user for download parameters
    
    Parameters:
    - options (argparse.Namespace): Command line arguments
    
    Returns:
    - dict: User input parameters
    """
    print("\nInteractive mode. Please provide necessary information.\n")

    # Database options
    database = options.database or validate_input(
        f"Enter database (genbank, refseq). Use commas for multiple entries [default: {DEFAULT_DATABASE}]: ",
        VALID_DATABASES, default=DEFAULT_DATABASE
    )

    # Organism options
    organism_group = options.organism_group or validate_input(
        "\nEnter organism group (archaea, bacteria, fungi, human, invertebrate, metagenomes, "
        "other, plant, protozoa, vertebrate_mammalian, vertebrate_other, viral). "
        "\nUse commas for multiple entries. Leave empty for all: ",
        VALID_ORGANISM_GROUPS, allow_empty=True
    )

    # Taxonomy ID
    taxid = options.taxid or input(
        "\nEnter taxonomy ID (e.g., 562 for NCBI or s__Escherichia coli for GTDB). "
        "Use commas for multiple entries. Leave empty for all: "
    )

    # File options
    file_types = options.file_types or validate_input(
        f"\nEnter file types (genomic.fna.gz, assembly_report.txt, protein.faa.gz, genomic.gbff.gz) "
        f"\nUse commas for multiple entries [default: {DEFAULT_FILE_TYPE}]: ",
        VALID_FILE_TYPES, default=DEFAULT_FILE_TYPE
    )

    # Filter options
    refseq_category = options.refseq_category or validate_input(
        f"\nEnter RefSeq category (reference genome, na). "
        f"\nUse commas for multiple entries [default: {DEFAULT_REFSEQ_CATEGORY}]: ",
        VALID_REFSEQ_CATEGORIES, default=DEFAULT_REFSEQ_CATEGORY
    )

    assembly_level = options.assembly_level or validate_input(
        f"\nEnter assembly level (complete genome, chromosome, scaffold, contig). "
        f"\nUse commas for multiple entries [default: {DEFAULT_ASSEMBLY_LEVEL}]: ",
        VALID_ASSEMBLY_LEVELS, default=DEFAULT_ASSEMBLY_LEVEL
    )

    start_date = options.start_date or input(
        "\nEnter start date (>=) for sequence release (format YYYYMMDD). Leave empty for no filter: "
    )

    end_date = options.end_date or input(
        "\nEnter end date (<=) for sequence release (format YYYYMMDD). Leave empty for no filter: "
    )

    custom_filter = options.custom_filter or input(
        "\nEnter custom filter for assembly summary (format colA:val1|colB:valX,valY). Leave empty for no filter: "
    )

    # Taxonomy options
    taxonomy_mode = options.taxonomy_mode or validate_input(
        f"\nEnter taxonomy mode (ncbi, gtdb) [default: {DEFAULT_TAXONOMY_MODE}]: ",
        VALID_TAXONOMY_MODES, default=DEFAULT_TAXONOMY_MODE
    )

    limit_assembly = options.limit_assembly or input(
        "\nEnter number of assemblies to download per taxa (0 for all, or rank:number e.g. genus:3) [default: 0]: ") or "0"

    # Run options
    output_dir = options.output_dir or input(
        f"\nEnter output directory [default: {DEFAULT_OUTPUT_DIR}]: ") or DEFAULT_OUTPUT_DIR

    threads = options.threads or input(f"\nEnter number of threads [default: {DEFAULT_THREADS}]: ") or DEFAULT_THREADS

    dry_run = options.dry_run or input("\nEnable dry run mode (no actual downloads) [y/N]: ").lower() == 'y'
    fix_mode = options.fix_mode or input("\nEnable fix mode (re-download incomplete or failed files) [y/N]: ").lower() == 'y'
    md5_check = options.md5_check or input("\nCheck MD5 of downloaded files [Y/n]: ").lower() != 'n'

    # Report options
    assembly_report = options.assembly_report or input("\nGenerate updated assembly accessions report [y/N]: ").lower() == 'y'
    sequence_report = options.sequence_report or input("\nGenerate updated sequence accessions report [y/N]: ").lower() == 'y'
    url_report = options.url_report or input("\nGenerate URLs report (success/failed) [y/N]: ").lower() == 'y'

    # Misc options
    version_label = options.version_label or input("\nEnter version label (leave empty for timestamp): ")
    external_assembly = options.external_assembly or input("\nEnter external assembly_summary.txt file path (leave empty for none): ")
    alt_version_label = options.alt_version_label or input("\nEnter alternative version label (leave empty for none): ")
    retry_attempts = options.retry_attempts or input(f"\nEnter number of retry attempts [default: {DEFAULT_RETRY_ATTEMPTS}]: ") or DEFAULT_RETRY_ATTEMPTS
    conditional_exit = options.conditional_exit or input("\nEnter conditional exit status (0 for off): ") or "0"
    ncbi_folders = options.ncbi_folders or input("\nOutput files in NCBI folder structure [y/N]: ").lower() == 'y'
    downloader = options.downloader or validate_input(
        f"\nEnter downloader to use (wget, curl) [default: {DEFAULT_DOWNLOADER}]: ",
        VALID_DOWNLOADERS, default=DEFAULT_DOWNLOADER
    )
    delete_extra = options.delete_extra or input("\nAllow deletion of extra files in output folder [y/N]: ").lower() == 'y'

    # Output mode
    silent = options.silent or input("\nEnable silent output [y/N]: ").lower() == 'y'
    progress_only = options.progress_only or input("\nEnable silent output with download progress only [y/N]: ").lower() == 'y'
    verbose = options.verbose or input("\nEnable verbose log [y/N]: ").lower() == 'y'
    debug = options.debug or input("\nEnable debug mode [y/N]: ").lower() == 'y'

    return {
        # Database options
        "database": database,
        # Organism options
        "organism_group": organism_group,
        "taxid": taxid,
        # File options
        "file_types": file_types,
        # Filter options
        "refseq_category": refseq_category,
        "assembly_level": assembly_level,
        "start_date": start_date,
        "end_date": end_date,
        "custom_filter": custom_filter,
        # Taxonomy options
        "taxonomy_mode": taxonomy_mode,
        "limit_assembly": limit_assembly,
        # Run options
        "output_dir": output_dir,
        "threads": threads,
        "dry_run": dry_run,
        "fix_mode": fix_mode,
        "md5_check": md5_check,
        # Report options
        "assembly_report": assembly_report,
        "sequence_report": sequence_report,
        "url_report": url_report,
        # Misc options
        "version_label": version_label,
        "external_assembly": external_assembly,
        "alt_version_label": alt_version_label,
        "retry_attempts": retry_attempts,
        "conditional_exit": conditional_exit,
        "ncbi_folders": ncbi_folders,
        "downloader": downloader,
        "delete_extra": delete_extra,
        # Output mode
        "silent": silent,
        "progress_only": progress_only,
        "verbose": verbose,
        "debug": debug
    }


def build_command(options):
    """
    Build genome_updater.sh command based on options
    
    Parameters:
    - options (argparse.Namespace): Command line arguments
    
    Returns:
    - str: Built command string
    """
    cmd = ["genome_updater.sh"]

    # Database options
    if options.database:
        cmd.append(f"-d '{options.database}'")

    # Organism options
    if options.organism_group:
        cmd.append(f"-g '{options.organism_group}'")

    if options.taxid:
        cmd.append(f"-T '{options.taxid}'")

    # File options
    if options.file_types:
        cmd.append(f"-f '{options.file_types}'")

    # Filter options
    if options.refseq_category:
        cmd.append(f"-c '{options.refseq_category}'")
        
    if options.assembly_level:
        cmd.append(f"-l '{options.assembly_level}'")
    
    if options.start_date:
        cmd.append(f"-D '{options.start_date}'")
    
    if options.end_date:
        cmd.append(f"-E '{options.end_date}'")
    
    if options.custom_filter:
        cmd.append(f"-F '{options.custom_filter}'")
    
    # Taxonomy options
    if options.taxonomy_mode:
        cmd.append(f"-M '{options.taxonomy_mode}'")
    
    if options.limit_assembly:
        # Check if limit_assembly contains a rank specification
        if ':' in str(options.limit_assembly):
            cmd.append(f"-A '{options.limit_assembly}'")
        else:
            cmd.append(f"-A 'species:{options.limit_assembly}'")
    
    # Keep taxonomy database
    cmd.append("-a")
    
    # Run options
    if options.output_dir:
        cmd.append(f"-o '{options.output_dir}'")

    if options.threads:
        cmd.append(f"-t {options.threads}")

    if options.dry_run:
        cmd.append("-k")

    if options.fix_mode:
        cmd.append("-i")

    if options.md5_check:
        cmd.append("-m")
    
    # Report options
    if getattr(options, 'assembly_report', False):
        cmd.append("-u")
    
    if getattr(options, 'sequence_report', False):
        cmd.append("-r")
    
    if getattr(options, 'url_report', False):
        cmd.append("-p")
    
    # Misc options
    if options.version_label:
        cmd.append(f"-b '{options.version_label}'")
    
    if options.external_assembly:
        cmd.append(f"-e '{options.external_assembly}'")
    
    if options.alt_version_label:
        cmd.append(f"-B '{options.alt_version_label}'")
    
    if options.retry_attempts:
        cmd.append(f"-R {options.retry_attempts}")
    
    if options.conditional_exit:
        cmd.append(f"-n {options.conditional_exit}")
    
    # NCBI folder structure
    if getattr(options, 'ncbi_folders', False):
        cmd.append("-N")
    
    if options.downloader:
        cmd.append(f"-L '{options.downloader}'")
    
    if getattr(options, 'delete_extra', False):
        cmd.append("-x")
    
    # Output mode settings
    if getattr(options, 'silent', False):
        cmd.append("-s")
    
    if getattr(options, 'progress_only', False):
        cmd.append("-w")
    
    if getattr(options, 'verbose', False):
        cmd.append("-V")
    
    if getattr(options, 'debug', False):
        cmd.append("-Z")

    return " ".join(cmd)


def run_command(cmd, ret_stdout=False, shell=False, quiet=False, debug=False, timeout=None):
    """
    Run command and handle results
    
    Parameters:
    - cmd (str or list): Command to execute
    - ret_stdout (bool): Whether to capture and return stdout
    - shell (bool): Whether to run command in shell
    - quiet (bool): Whether to suppress stderr output
    - debug (bool): Whether to enable debug mode
    - timeout (int): Command execution timeout in seconds
    
    Returns:
    - str: Captured stdout if ret_stdout is True, otherwise None
    """
    errcode = 1
    stdout = None

    try:
        if debug:
            print(f"Executing command: {cmd}")

        # Split command string into list if not using shell
        if isinstance(cmd, str) and not shell:
            cmd = shlex.split(cmd)

        # Run command
        process = subprocess.Popen(
            cmd,
            shell=shell,
            universal_newlines=True,
            stdout=subprocess.PIPE if ret_stdout else sys.stderr,
            stderr=subprocess.PIPE if quiet else sys.stderr
        )

        if ret_stdout:
            stdout, stderr = process.communicate(timeout=timeout)
            if stderr and not quiet:
                print(stderr)
        else:
            process.wait(timeout=timeout)

        errcode = process.returncode
        if errcode != 0:
            raise subprocess.CalledProcessError(errcode, cmd)

    except subprocess.CalledProcessError as e:
        print(f"Command execution failed with code {e.returncode}: {e.cmd}")
        sys.exit(e.returncode)
    except subprocess.TimeoutExpired:
        print(f"Command timed out after {timeout} seconds: {cmd}")
        process.kill()
        sys.exit(1)
    except Exception as e:
        print(f"Unexpected error running command: {cmd}")
        print(str(e))
        sys.exit(errcode)

    return stdout


def validate_input_files(input_files_folder, file_types, quiet=False, input_recursive=False):
    """
    Validate input files
    
    Parameters:
    - input_files_folder (list): List of input files or folders
    - file_types (list): List of file types
    - quiet (bool): Whether to suppress output
    - input_recursive (bool): Whether to search folders recursively
    
    Returns:
    - set: Set of valid files
    """
    valid_input_files = set()

    for i in input_files_folder:
        if os.path.isfile(i):
            valid_input_files.add(i)
        elif os.path.isdir(i):
            files_in_dir = 0

            if input_recursive:
                for file_type in file_types:
                    for path in Path(i).rglob('*' + file_type):
                        f = str(path)
                        if os.path.isfile(f):
                            files_in_dir += 1
                            valid_input_files.add(f)
            else:
                for file in os.listdir(i):
                    for file_type in file_types:
                        if file.endswith(file_type):
                            f = os.path.join(i, file)
                            if os.path.isfile(f):
                                files_in_dir += 1
                                valid_input_files.add(f)

            if not quiet:
                print(f" - {files_in_dir} valid files [types: {', '.join(file_types)}" +
                      (", recursive" if input_recursive else "") +
                      f"] found in {i}")
        else:
            if not quiet:
                print(f" - Skipping invalid file/folder: {i}")

    if not quiet:
        print(f" - Total valid files: {len(valid_input_files)}")

    return valid_input_files


def load_assembly_accession(input_files):
    """
    Load assembly accessions from input files
    
    Parameters:
    - input_files (set): Set of input files
    
    Returns:
    - pandas.DataFrame: DataFrame containing assembly accession information
    """
    info_cols = ["file", "target", "node", "specialization", "specialization_name"]

    # Regular expression to match GenBank/RefSeq assembly accessions
    pattern = re.compile(r"GC[AF]_[0-9]+\.[0-9]+")

    # Extract assembly accession or use filename as target
    data = []
    for f in input_files:
        match = pattern.search(f)
        target = match.group() if match else os.path.basename(f)
        data.append((target, f))

    # Create DataFrame with info_cols columns
    info = pd.DataFrame(data, columns=["target", "file"])

    # Add other columns with default values
    for col in info_cols:
        if col not in info.columns:
            info[col] = None

    # Clean data
    info.dropna(how="all", inplace=True)
    info.dropna(subset=["target"], inplace=True)
    info.drop_duplicates(subset=["target"], inplace=True)

    # Set 'target' as index
    info.set_index('target', inplace=True)

    return info


def parse_assembly_summary(info, assembly_summary, level="species"):
    """
    Parse assembly summary file
    
    Parameters:
    - info (pandas.DataFrame): Assembly information DataFrame
    - assembly_summary (str): Path to assembly summary file
    - level (str): Taxonomy level
    
    Returns:
    - dict: Assembly summary counts
    """
    count_assembly_summary = {}
    unique_acc = set(info.index)

    # Detect header lines
    header_lines = 0
    with open(assembly_summary, 'r') as ass_sum:
        for line in ass_sum:
            if line[0] == "#":
                header_lines += 1
            else:
                break

    # Read assembly summary file
    tmp_acc_node = pd.read_csv(
        assembly_summary,
        sep="\t",
        header=None,
        skiprows=header_lines,
        usecols=[0, 5, 7, 8],
        names=["target", "node", "organism_name", "infraspecific_name"],
        index_col="target",
        converters={"target": lambda x: x if x in unique_acc else None, "node": str}
    )
    
    # Keep only used sequence IDs
    tmp_acc_node = tmp_acc_node[tmp_acc_node.index.notnull()]

    # Save count
    count_assembly_summary[assembly_summary] = tmp_acc_node.shape[0]

    # Create specialization
    if level == "assembly":
        # Handle infraspecific_name prefix
        tmp_acc_node["infraspecific_name"] = tmp_acc_node["infraspecific_name"].replace(
            "^[a-z]+=", "", regex=True).fillna("")

        # Build name
        def build_name(n):
            if n.organism_name.endswith(n.infraspecific_name):
                return n.organism_name
            else:
                return n.organism_name + " " + n.infraspecific_name

        # Add infraspecific_name suffix
        tmp_acc_node["specialization_name"] = tmp_acc_node[["organism_name", "infraspecific_name"]].apply(
            lambda n: build_name(n), axis=1)
        tmp_acc_node["specialization"] = tmp_acc_node.index

    # Merge nodes and specializations found by target (accession)
    if count_assembly_summary[assembly_summary]:
        info.update(tmp_acc_node)
    
    del tmp_acc_node
    return count_assembly_summary


def get_file_info(options, info, build_output_folder, assembly_summary):
    """
    Get file information
    
    Parameters:
    - options (argparse.Namespace): Command line arguments
    - info (pandas.DataFrame): Assembly information DataFrame
    - build_output_folder (str): Build output folder
    - assembly_summary (str): Path to assembly summary file
    """
    # Parse assembly summary file
    start_time = time.time()
    quiet = getattr(options, 'silent', False) or getattr(options, 'progress_only', False)
    
    if not quiet:
        print("Parsing assembly summary file")
    
    count_assembly_summary = parse_assembly_summary(info, assembly_summary)

    # Output parsing results
    if not quiet:
        for assembly_summary_file, count in count_assembly_summary.items():
            file_name = assembly_summary_file.split("/")[-1]
            print(f" - Found {count} entries in {file_name}")

        print(f" - Done in {time.time() - start_time:.2f} seconds.\n")


def validate_taxonomy(info, tax):
    """
    Validate taxonomy: convert to species level nodes
    
    Parameters:
    - info (pandas.DataFrame): Assembly information DataFrame
    - tax (multitax.NcbiTx): Taxonomy object
    """
    # Convert nodes to species level parent nodes
    info["node"] = info["node"].apply(lambda n: tax.parent_rank(n, "species"))

    # Skip invalid nodes
    na_entries = info["node"].isna().sum()
    if na_entries > 0:
        info.dropna(subset=["node"], inplace=True)


def remove_files(files):
    """
    Remove files
    
    Parameters:
    - files (str or list): Files to remove
    """
    if isinstance(files, str):
        files = [files]
    
    for f in files:
        if os.path.isfile(f):
            os.remove(f)


def write_tax(tax_file, info, tax):
    """
    Write taxonomy file
    
    Parameters:
    - tax_file (str): Path to taxonomy file
    - info (pandas.DataFrame): Assembly information DataFrame
    - tax (multitax.NcbiTx): Taxonomy object
    """
    # Remove existing file
    remove_files(tax_file)
    
    # Write taxonomy to file
    tax.write(tax_file)

    # Load written taxonomy file
    tax_df = pd.read_csv(
        tax_file, 
        names=["node", "parent", "rank", "name"], 
        delimiter='\t', 
        dtype=str
    )

    # Write back updated DataFrame
    tax_df.to_csv(tax_file, sep="\t", header=False, index=False)

    # Create target.tsv file
    target_tsv_path = os.path.join(os.path.dirname(tax_file), "target.tsv")

    # Keep only file and node columns
    info_subset = info[["file", "node"]]

    # Write info_subset content to target_tsv_path
    with open(target_tsv_path, 'w') as f:
        f.write(info_subset.to_csv(sep="\t", header=False, index=False))


def setup_output_directory(options):
    """
    Set up output directory
    
    Parameters:
    - options (argparse.Namespace): Command line arguments
    
    Returns:
    - tuple: (output_folder, assembly_summary, tmp_folder)
    """
    output_folder = options.output_dir
    
    # Handle output directory
    if options.fix_mode:
        print("Fix mode enabled. Will re-download existing files.")
    elif os.path.exists(output_folder):
        while True:
            confirmation = input(
                f"Output folder '{output_folder}' exists. Do you want to clear it [y], continue [c], or cancel [n]? [y/c/n]: ").lower()
            
            if confirmation == 'y':
                print("Clearing output folder...")
                if os.name == 'nt':  # Windows
                    os.system(f"rmdir /s /q {output_folder}")
                else:  # Unix/Linux
                    os.system(f"rm -rf {output_folder}")
                os.makedirs(output_folder)  # Recreate folder after clearing
                break
            elif confirmation == 'c':
                print("Continuing without clearing output folder.")
                break
            elif confirmation == 'n':
                print("Operation cancelled.")
                sys.exit(1)
            else:
                print("Invalid option. Please choose 'y', 'c', or 'n'.")
    else:
        print(f"Creating output directory: {output_folder}")
        os.makedirs(output_folder, exist_ok=True)

    # Set up paths
    assembly_summary = os.path.join(output_folder, "assembly_summary.txt")
    tmp_folder = os.path.join(output_folder, "tmp")
    os.makedirs(tmp_folder, exist_ok=True)

    # Check existing files
    if os.path.isfile(assembly_summary) and os.path.getsize(assembly_summary) > 0:
        print(f"Assembly summary file exists: {assembly_summary}")

    return output_folder, assembly_summary, tmp_folder


def download(interactive=False, raw_args=None):
    """
    Download NCBI genome data
    
    Parameters:
    - interactive (bool): Whether to use interactive mode
    
    Returns:
    - argparse.Namespace: Command line arguments
    """
    # Create argument parser
    parser = argparse.ArgumentParser(
        description='Download NCBI genome data',
        epilog='''Examples:
        python download.py -d refseq -g "archaea,bacteria" -l "complete genome" -f genomic.fna.gz -o ./output
        ''',
        formatter_class=argparse.RawTextHelpFormatter
    )

    # Add arguments
    # Database options
    parser.add_argument("-d", "--database", 
                       help=f"Database to download (default: {DEFAULT_DATABASE})\nOptions: {', '.join(VALID_DATABASES)}")
    
    # Organism options
    parser.add_argument("-g", "--organism-group",
                       help=f"Organism groups to download (comma-separated)\nOptions: {', '.join(VALID_ORGANISM_GROUPS)}")
    
    parser.add_argument("-T", "--taxid",
                       help="Taxonomy IDs to download (comma-separated)\nExample: 562,623 (for -M ncbi) or s__Escherichia coli (for -M gtdb)")
    
    # File options
    parser.add_argument("-f", "--file-types",
                       help=f"File types to download (default: {DEFAULT_FILE_TYPE})\nOptions: {', '.join(VALID_FILE_TYPES)}")
    
    # Filter options
    parser.add_argument("-c", "--refseq-category",
                       help=f"RefSeq categories to download (default: {DEFAULT_REFSEQ_CATEGORY})\nOptions: {', '.join(VALID_REFSEQ_CATEGORIES)}")
    
    parser.add_argument("-l", "--assembly-level",
                       help=f"Assembly levels to download (default: {DEFAULT_ASSEMBLY_LEVEL})\nOptions: {', '.join(VALID_ASSEMBLY_LEVELS)}")
    
    parser.add_argument("-D", "--start-date",
                       help="Start date (>=), based on sequence release date. Format YYYYMMDD")
    
    parser.add_argument("-E", "--end-date",
                       help="End date (<=), based on sequence release date. Format YYYYMMDD")
    
    parser.add_argument("-F", "--custom-filter",
                       help="Custom filter for assembly summary in format colA:val1|colB:valX,valY (case insensitive)")
    
    # Taxonomy options
    parser.add_argument("-M", "--taxonomy-mode",
                       help=f"Taxonomy mode (default: {DEFAULT_TAXONOMY_MODE})\nOptions: {', '.join(VALID_TAXONOMY_MODES)}")
    
    parser.add_argument("-A", "--limit-assembly",
                       help="Limit number of assemblies for each selected taxa. 0 for all.\nSelection by ranks also supported with rank:number (e.g., genus:3)")
    
    # Run options
    parser.add_argument("-o", "--output-dir",
                       help=f"Output directory path (default: {DEFAULT_OUTPUT_DIR})\nWarning: Will be cleared if exists")
    
    parser.add_argument("-t", "--threads",
                       help=f"Number of download threads (default: {DEFAULT_THREADS})")
    
    parser.add_argument("-k", "--dry-run",
                       action="store_true",
                       help="Enable dry run mode (no actual downloads)")
    
    parser.add_argument("-i", "--fix-mode",
                       action="store_true", 
                       help="Only re-download incomplete or failed files")
    
    parser.add_argument("-m", "--md5-check",
                       action="store_true",
                       help="Check MD5 of downloaded files")
    
    # Report options
    parser.add_argument("-u", "--assembly-report",
                       action="store_true",
                       help="Generate updated assembly accessions report")
    
    parser.add_argument("-r", "--sequence-report",
                       action="store_true",
                       help="Generate updated sequence accessions report")
    
    parser.add_argument("-p", "--url-report",
                       action="store_true",
                       help="Generate URLs report (success/failed)")
    
    # Misc options
    parser.add_argument("-b", "--version-label",
                       help="Version label (default: current timestamp)")
    
    parser.add_argument("-e", "--external-assembly",
                       help="External assembly_summary.txt file to recover data from")
    
    parser.add_argument("-B", "--alt-version-label",
                       help="Alternative version label to use as current version")
    
    parser.add_argument("-R", "--retry-attempts",
                       help=f"Number of attempts to retry download files in batches (default: {DEFAULT_RETRY_ATTEMPTS})")

    parser.add_argument("-n", "--conditional-exit",
                        help="Conditional exit status based on number of failures accepted")

    parser.add_argument("-N", "--ncbi-folders",
                        action="store_true",
                        help="Output files in folders like NCBI ftp structure")

    parser.add_argument("-L", "--downloader",
                        help=f"Downloader to use (default: {DEFAULT_DOWNLOADER})\nOptions: {', '.join(VALID_DOWNLOADERS)}")

    parser.add_argument("-x", "--delete-extra",
                        action="store_true",
                        help="Allow deletion of regular extra files found in output folder")

    parser.add_argument("-s", "--silent",
                        action="store_true",
                        help="Silent output")

    parser.add_argument("-w", "--progress-only",
                        action="store_true",
                        help="Silent output with download progress only")

    parser.add_argument("-V", "--verbose",
                        action="store_true",
                        help="Verbose log")

    parser.add_argument("-Z", "--debug",
                        action="store_true",
                        help="Print debug information and run in debug mode")

    
    if raw_args is not None:
        options = parser.parse_args(raw_args)
        interactive = False
    elif interactive:
        parsed_args = parser.parse_args([])
        options = argparse.Namespace(**prompt_user(parsed_args))
    else:
        # 命令行模式
        options = parser.parse_args()
    
    # Build command
    command = build_command(options)

    # Set up output directory
    output_folder, assembly_summary, tmp_folder = setup_output_directory(options)

    # Execute download
    try:
        # Download phase
        quiet = getattr(options, 'silent', False) or getattr(options, 'progress_only', False)
        if not quiet:
            download_start = time.time()
            print("Downloading data...")
        
        # Execute download command
        run_command(command, shell=True, quiet=quiet)
        
        if not quiet:
            download_end = time.time()
            print(f"\nDownload completed in {download_end - download_start:.2f} seconds.")

        # Processing phase
        process_start = time.time()
        if not quiet:
            print("\nProcessing data...")
        
        # Get file paths
        file_folder = os.path.join(options.output_dir, os.path.dirname(os.readlink(assembly_summary)), "files")
        input_file = validate_input_files([file_folder], options.file_types.split(","), quiet, input_recursive=True)

        if not input_file:
            if not quiet:
                print("No valid files found. Exiting.")
            sys.exit(1)
        else:
            if not quiet:
                print("Valid files found.")

        # Load taxonomy information
        tax = NcbiTx()
        info = load_assembly_accession(input_file)
        
        if info.empty:
            if not quiet:
                print("No valid assembly accessions found. Exiting.")
            sys.exit(1)
        
        # Get file information and validate taxonomy
        get_file_info(options, info, tmp_folder, assembly_summary)
        validate_taxonomy(info, tax)
        
        # Filter taxonomy
        unique_nodes = info["node"].unique()
        tax.filter(unique_nodes)
        
        # Write taxonomy file
        tax_file_path = os.path.join(output_folder, "tax.info")
        write_tax(tax_file_path, info, tax)
        
        if not quiet:
            process_end = time.time()
            print(f"\nProcessing completed in {process_end - process_start:.2f} seconds.")

        # Clean up temporary folder
        if os.name == 'nt':  # Windows
            os.system(f"rmdir /s /q {tmp_folder}")
        else:  # Unix/Linux
            os.system(f"rm -rf {tmp_folder}")

    except Exception as e:
        print(f"Error during download process: {str(e)}")
        sys.exit(1)

    return options

def main():
    """
    Main function to execute when script is run directly
    """
    try:
        # Auto-detect mode: use argument mode if command line arguments exist, otherwise use interactive mode
        options = download(interactive=None)
        print("\nDownload completed!")
        print(f"Output directory: {options.output_dir}")
    except KeyboardInterrupt:
        print("\nOperation interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\nError occurred: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()