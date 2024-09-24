import argparse
import re
import subprocess
import sys
import time
import os
import shlex
import urllib
from pathlib import Path
import pandas as pd
from multitax import NcbiTx


def prompt_user(options):
    print("\nInteractive mode. Please provide the necessary information.\n")

    valid_databases = ["genbank", "refseq"]
    valid_organism_groups = ["archaea", "bacteria", "fungi", "human", "invertebrate", "metagenomes",
                             "other", "plant", "protozoa", "vertebrate_mammalian", "vertebrate_other", "viral"]
    valid_assembly_levels = ["complete genome", "chromosome", "scaffold", "contig"]
    valid_refseq_categories = ["reference genome", "representative genome", "na"]
    valid_file_types = ["genomic.fna.gz", "assembly_report.txt", "protein.faa.gz", "genomic.gbff.gz"]

    def validate_input(prompt, valid_options, default=None, allow_empty=False):
        while True:
            user_input = input(prompt).strip()
            if not user_input and default is not None:
                return default
            if allow_empty and not user_input:
                return ""
            entries = [entry.strip() for entry in user_input.split(",")]
            if all(entry in valid_options for entry in entries):
                return ",".join(entries)
            print(f"\nInvalid input. Please enter a valid option from: {', '.join(valid_options)}")

    database = options.database or validate_input(
        "Enter the database(s) (genbank, refseq). Use comma to separate multiple entries [default: refseq]: ",
        valid_databases, default="refseq"
    )

    organism_group = options.organism_group or validate_input(
        "\nEnter organism group(s) (archaea, bacteria, fungi, human, invertebrate, metagenomes, "
        "other, plant, protozoa, vertebrate_mammalian, vertebrate_other, viral). "
        "\nUse comma to separate multiple entries. Leave empty for all: ",
        valid_organism_groups, allow_empty=True
    )

    taxid = options.taxid or input(
        "\nEnter taxonomy ID(s) (e.g., 562). Use comma to separate multiple entries. Leave empty for all: "
    )

    assembly_level = options.assembly_level or validate_input(
        "\nEnter assembly level(s) (complete genome, chromosome, scaffold, contig). "
        "\nUse comma to separate multiple entries. [default: complete genome]: ",
        valid_assembly_levels, default="complete genome"
    )

    refseq_category = options.refseq_category or validate_input(
        "\nEnter RefSeq category (reference genome, representative genome, na). "
        "\nUse comma to separate multiple entries [default: representative genome]: ",
        valid_refseq_categories, default="representative genome"
    )

    file_types = options.file_types or validate_input(
        "\nEnter file type(s) (genomic.fna.gz, assembly_report.txt, protein.faa.gz, genomic.gbff.gz) "
        "\nUse comma to separate multiple entries [default: genomic.fna.gz]: ",
        valid_file_types, default="genomic.fna.gz"
    )

    limit_assembly = options.limit_assembly or int(
        input("\nEnter the number of assemblies to download [default: 0]: ") or 0)

    output_dir = options.output_dir or input(
        "\nEnter the output directory [default: ./genome_output]: ") or "./genome_output"

    threads = options.threads or input("\nEnter the number of threads [default: 1]: ") or "1"

    dry_run = options.dry_run or input("\nEnable dry-run mode (no data will be downloaded) [y/N]: ").lower() == 'y'

    fix_mode = options.fix_mode or input(
        "\nEnable fix-only mode (re-download incomplete or failed data) [y/N]: ").lower() == 'y'

    md5_check = options.md5_check or input("\nCheck MD5 of downloaded files [Y/n]: ").lower() != 'n'

    quite_or_verbose = input("\nquite or verbose mode [q/V]: ").lower() == 'q'

    return {
        "database": database,
        "organism_group": organism_group,
        "file_types": file_types,
        "taxid": taxid,
        "assembly_level": assembly_level,
        "refseq_category": refseq_category,
        "limit_assembly": limit_assembly,
        "output_dir": output_dir,
        "threads": threads,
        "dry_run": dry_run,
        "fix_mode": fix_mode,
        "md5_check": md5_check,
        "quite_or_verbose": quite_or_verbose
    }


def build_command(options):
    cmd = ["genome_updater.sh"]

    if options.database:
        cmd.append(f"-d '{options.database}'")

    if options.organism_group:
        cmd.append(f"-g '{options.organism_group}'")

    if options.taxid:
        cmd.append(f"-T '{options.taxid}'")

    if options.file_types:
        cmd.append(f"-f '{options.file_types}'")

    if options.assembly_level:
        cmd.append(f"-l '{options.assembly_level}'")

    if options.refseq_category:
        cmd.append(f"-c '{options.refseq_category}'")

    if options.limit_assembly:
        cmd.append(f"-A 'species:{options.limit_assembly}'")

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

    if options.quite_or_verbose:
        cmd.append("-q -w")
    else:
        cmd.append("-V")

    cmd.append("-N")

    return " ".join(cmd)


def run(cmd, ret_stdout: bool = False, shell: bool = False, quiet: bool = False, debug: bool = False,
        timeout: int = None):
    """
    Runs a command with optional configurations.

    Parameters:
    - cmd (str or list): Command to execute. Can be a string or a list.
    - ret_stdout (bool): Capture and return stdout if True.
    - shell (bool): Run the command in a shell if True.
    - quiet (bool): Suppress stderr if True.
    - debug (bool): Enable debug mode for more detailed output.
    - timeout (int): Timeout for the command execution in seconds.

    Returns:
    - stdout (str): Captured stdout if ret_stdout is True, otherwise None.
    """
    errcode = 1
    stdout = None

    try:
        if debug:
            print(f"Executing command: {cmd}")

        # If cmd is a string and not using shell, split it into a list
        if isinstance(cmd, str) and not shell:
            cmd = shlex.split(cmd)

        # Run the command
        process = subprocess.Popen(cmd if shell else cmd,
                                   shell=shell,
                                   universal_newlines=True,  # text= (from py3.7)
                                   stdout=subprocess.PIPE if ret_stdout else sys.stderr,
                                   stderr=subprocess.PIPE if quiet else sys.stderr)

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
        print(f"The command failed with error code {e.returncode}: {e.cmd}")
        sys.exit(e.returncode)
    except subprocess.TimeoutExpired:
        print(f"The command timed out after {timeout} seconds: {cmd}")
        process.kill()
        sys.exit(1)
    except Exception as e:
        print(f"An unexpected error occurred while running the command: {cmd}")
        print(str(e))
        sys.exit(errcode)

    return stdout


def validate_input_files(input_files_folder, file_types, quiet, input_recursive: bool = False):
    """
    Given a list of input files and/or folders and a list of file types,
    check for valid files and return them in a set.
    """
    valid_input_files = set()  # 用于存储有效文件的集合

    for i in input_files_folder:  # 遍历输入的文件或文件夹
        if os.path.isfile(i):  # 检查 i 是否为有效文件
            valid_input_files.add(i)  # 如果是有效文件，添加到集合中
        elif os.path.isdir(i):  # 如果 i 是文件夹
            files_in_dir = 0  # 记录找到的有效文件数量

            if input_recursive:  # 如果设置为递归查找
                for file_type in file_types:  # 遍历所有文件类型
                    for path in Path(i).rglob('*' + file_type):  # 递归查找符合文件类型的文件
                        f = str(path)
                        if os.path.isfile(f):  # 检查文件是否存在
                            files_in_dir += 1
                            valid_input_files.add(f)  # 添加到集合中
            else:  # 非递归查找，只检查顶层目录
                for file in os.listdir(i):  # 遍历文件夹中的文件
                    for file_type in file_types:  # 遍历所有文件类型
                        if file.endswith(file_type):  # 过滤符合文件类型的文件
                            f = os.path.join(i, file)
                            if os.path.isfile(f):  # 检查文件是否存在
                                files_in_dir += 1
                                valid_input_files.add(f)  # 添加到集合中

            # 打印找到的文件数量和其他相关信息
            if not quiet:
                print(f" - {files_in_dir} valid file(s) [file types: {', '.join(file_types)}" +
                      (", recursive" if input_recursive else "") +
                      f"] found in {i}")
        else:  # 如果 i 既不是有效文件也不是文件夹
            if not quiet:
                print(f" - skipping invalid file/folder: {i}")

    # 打印总的有效文件数量
    if not quiet:
        print(f" - total valid files: {len(valid_input_files)}")

    return valid_input_files  # 返回有效文件集合


def load_assembly_accession(input_files):
    info_cols = ["file", "target", "node", "specialization", "specialization_name"]

    # 正则表达式匹配 GenBank/RefSeq 的组装编号
    pattern = re.compile(r"GC[AF]_[0-9]+\.[0-9]+")

    # 提取组装编号或使用文件名作为目标
    data = [(pattern.search(f).group() if pattern.search(f) else os.path.basename(f), f) for f in input_files]

    # 创建包含 info_cols 列的数据框
    info = pd.DataFrame(data, columns=["target", "file"])

    # 添加缺省值的其他列
    for col in info_cols:
        if col not in info.columns:
            info[col] = None

    # 清理数据
    info.dropna(how="all", inplace=True)
    info.dropna(subset=["target"], inplace=True)
    info.drop_duplicates(subset=["target"], inplace=True)

    # 设置 'target' 为索引
    info.set_index('target', inplace=True)

    return info


def parse_assembly_summary(info, assembly_summary, level):
    count_assembly_summary = {}
    unique_acc = set(info.index)

    # Detect header lines beforehand, so pandas.read_csv can read it properly
    header_lines = 0
    with open(assembly_summary, 'r') as ass_sum:
        for line in ass_sum:
            if line[0] == "#":
                header_lines += 1
            else:
                break

    tmp_acc_node = pd.read_csv(assembly_summary,
                               sep="\t",
                               header=None,
                               skiprows=header_lines,
                               # usecols = 1:assembly_accession, 6:taxid, 8:organism_name, 9:infraspecific_name
                               usecols=[0, 5, 7, 8],
                               names=["target", "node", "organism_name", "infraspecific_name"],
                               index_col="target",
                               converters={"target": lambda x: x if x in unique_acc else None, "node": str})
    tmp_acc_node = tmp_acc_node[tmp_acc_node.index.notnull()]  # keep only seqids used

    # save count to return
    count_assembly_summary[assembly_summary] = tmp_acc_node.shape[0]

    # Create specialization
    if level == "assembly":
        # infraspecific_name has a prefix: breed=, cultivar=, ecotype= or strain=
        tmp_acc_node["infraspecific_name"] = tmp_acc_node["infraspecific_name"].replace("^[a-z]+=", "",
                                                                                        regex=True).fillna("")

        def build_name(n):
            if n.organism_name.endswith(n.infraspecific_name):
                return n.organism_name
            else:
                return n.organism_name + " " + n.infraspecific_name

        # add sufix of the infraspecific_name if not yet contained in the end of organism_name
        tmp_acc_node["specialization_name"] = tmp_acc_node[["organism_name",
                                                            "infraspecific_name"]].apply(lambda n: build_name(n),
                                                                                         axis=1)
        tmp_acc_node["specialization"] = tmp_acc_node.index

    # merge node(taxid) and specialization retrieved based on target(accesion)
    if count_assembly_summary[assembly_summary]:
        info.update(tmp_acc_node)
    del tmp_acc_node

    return count_assembly_summary


def get_file_info(options, info, build_output_folder, assemblySummary):
    """
    Load and parse NCBI assembly_summary files based on the provided options.
    """

    # 解析 assembly_summary 文件
    start_time = time.time()
    if not options.quite_or_verbose:
        print("Parsing assembly_summary files")
    count_assembly_summary = parse_assembly_summary(info, assemblySummary, "species")

    # 输出解析结果
    if not options.quite_or_verbose:
        for assembly_summary_file, count in count_assembly_summary.items():
            file_name = assembly_summary_file.split("/")[-1]
            print(f" - {count} entries found in the {file_name} file")

        print(f" - done in {time.time() - start_time:.2f}s.\n")


def validate_taxonomy(info, tax, options):
    """
    Validate taxonomy: convert to species level nodes.
    """

    # 将节点转换为物种级别的父节点
    info["node"] = info["node"].apply(lambda n: tax.parent_rank(n, "species"))

    # 跳过无效节点 (na == tax.undefined_node (None))
    na_entries = info["node"].isna().sum()

    if na_entries > 0:
        info.dropna(subset=["node"], inplace=True)


def rm_files(files):
    if isinstance(files, str):
        files = [files]
    for f in files:
        if os.path.isfile(f):
            os.remove(f)


def write_tax(tax_file, info, tax):
    """
    Write tabular taxonomy file (.tax)
    May include specialization as nodes.
    """
    tax_rank = "species"
    for target, row in info.iterrows():
        tax_node = target
        tax_name = target

        # # Check if the node is already present with the correct parent
        # if tax.latest(tax_node) is tax.undefined_node:
        #     tax.add(tax_node, row["node"], name=tax_name, rank=tax_rank)
        # else:
        #     assert tax.parent(tax_node) == row["node"]

    # Write filtered taxonomy with added nodes
    rm_files(tax_file)  # Remove the existing file if present
    tax.write(tax_file)  # Write the taxonomy to file

    # Load the written taxonomy file into a DataFrame for further processing
    tax_df = pd.read_csv(tax_file, names=["node", "parent", "rank", "name"], delimiter='\t', dtype=str)

    # Write the updated DataFrame back to the file
    tax_df.to_csv(tax_file, sep="\t", header=False, index=False)

    target_tsv_path = os.path.join(os.path.dirname(tax_file), "target.tsv")

    # 只保留 file 和 node 列
    info_subset = info[["file", "node"]]

    # 把 info_subset 的内容写入 target_tsv_path，不保留列名
    with open(target_tsv_path, 'w') as f:
        f.write(info_subset.to_csv(sep="\t", header=False, index=False))


def download(interactive=False):
    parser = argparse.ArgumentParser(description='Download NCBI fasta file ')

    parser.add_argument("-d", "--database", help="Database (genbank, refseq)")
    parser.add_argument("-g", "--organism-group", help="Organism group(s) (e.g., bacteria, fungi, human)")
    parser.add_argument("-T", "--taxid", help="Taxonomy ID(s) (e.g., 562)")
    parser.add_argument("-f", "--file-types", help="File types (e.g., genomic.fna.gz, assembly_report.txt)")
    parser.add_argument("-l", "--assembly-level", help="Assembly level (e.g., complete genome, chromosome)")
    parser.add_argument("-c", "--refseq-category",
                        help="Refseq category (e.g., reference genome, representative genome, na)")
    parser.add_argument("-A", "--limit-assembly", help="Limit the number of assemblies to download")
    parser.add_argument("-o", "--output-dir",
                        help="Name of the output directory (required). WARNING: If the directory exists, it will be cleared.")
    parser.add_argument("-t", "--threads", help="Number of threads")
    parser.add_argument("-k", "--dry-run", action="store_true", help="Enable dry-run mode (no data will be downloaded)")
    parser.add_argument("-i", "--fix-mode", action="store_true",
                        help="Enable fix-only mode (re-download incomplete or failed data)")
    parser.add_argument("-m", "--md5-check", action="store_true", help="Check MD5 of downloaded files")

    parser.add_argument('-q', '--quit', action='store_true', help="Enable quit mode")
    parser.add_argument('-v', '--verbose', action='store_true', help="Enable verbose mode")

    if not interactive:
        args = parser.parse_args()
    else:
        args = parser.parse_args([])

    if not any(vars(args).values() or interactive):
        options = argparse.Namespace(**prompt_user(args))
    else:
        options = args
        if not options.quit:
            options.quite_or_verbose = False
    command = build_command(options)

    outputFolder = options.output_dir
    if options.fix_mode:
        print("Fix mode enabled. Existing files will be re-downloaded.")
    elif os.path.exists(outputFolder):
        while True:
            confirmation = input(
                f"Output folder '{outputFolder}' already exists. Do you want to clear it [y], continue [c], or cancel [n]? [y/c/n]: ").lower()
            if confirmation == 'y':
                print("Clearing the output folder...")
                os.system(f"rm -rf {outputFolder}")
                os.makedirs(outputFolder)  # 清除后重新创建文件夹
                break
            elif confirmation == 'c':
                print("Continuing without clearing the output folder.")
                break
            elif confirmation == 'n':
                print("Operation cancelled.")
                sys.exit(1)
            else:
                print("Invalid option. Please choose 'y', 'c', or 'n'.")
    else:
        print(f"Creating output directory: {outputFolder}")
        os.makedirs(outputFolder, exist_ok=True)

    assemblySummary = outputFolder + "/assembly_summary.txt"
    tmpFolder = outputFolder + "/tmp"
    os.makedirs(tmpFolder, exist_ok=True)

    if os.path.isfile(assemblySummary) and os.path.getsize(assemblySummary) > 0:
        print(f"Assembly summary file already exists: {assemblySummary}")

    if not options.quite_or_verbose:
        downloadStart = time.time()
        print("Downloading data...")
    run(command, shell=True, quiet=options.quite_or_verbose, ret_stdout=options.quite_or_verbose)
    if not options.quite_or_verbose:
        downloadEnd = time.time()
        print("\nDownload time: " + str("%.2f" % (downloadEnd - downloadStart)) + " seconds.")

    processStart = time.time()
    if not options.quite_or_verbose:
        print("\nProcessing data...")
    file_folder = os.path.join(options.output_dir, os.path.dirname(os.readlink(assemblySummary)), "files")
    input_file = validate_input_files([file_folder], options.file_types.split(","), options.quite_or_verbose,
                                      input_recursive=True)

    if input_file:
        if not options.quite_or_verbose:
            print("Valid files found.")
    else:
        if not options.quite_or_verbose:
            print("No valid files found. Exiting.")
        sys.exit(1)
    tax = NcbiTx()
    info = load_assembly_accession(input_file)
    if info.empty:
        if not options.quite_or_verbose:
            print("No valid assembly accessions found. Exiting.")
        sys.exit(1)
    get_file_info(options, info, tmpFolder, assemblySummary)
    validate_taxonomy(info, tax, options)
    unique_nodes = info["node"].unique()
    tax.filter(unique_nodes)
    tax_file_path = os.path.join(outputFolder, "tax.info")
    write_tax(tax_file_path, info, tax)
    processEnd = time.time()
    if not options.quite_or_verbose:
        print("\nProcessing time: " + str("%.2f" % (processEnd - processStart)) + " seconds.")

    os.system(f"rm -rf {tmpFolder}")

    return options


if __name__ == "__main__":
    download()
