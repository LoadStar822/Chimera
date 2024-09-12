# Chimera: A Versatile Metagenomic Classification Tool

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE) 
[![Build Status](https://img.shields.io/badge/build-passing-brightgreen)](https://github.com/MalabZ/Chimera)
[![Conda](https://img.shields.io/conda/vn/MALAB/chimera)](https://anaconda.org/MALAB/chimera)
[![Platform](https://img.shields.io/badge/platform-linux--64-lightgrey)](https://www.linux.org/)
[![Docker Pulls](https://img.shields.io/docker/pulls/tianqinzhong/chimera)](https://hub.docker.com/r/tianqinzhong/chimera)


## Table of Contents
- [Project Overview](#project-overview)
- [Installation](#installation)
   1. [Source Installation](#1-source-installation)
   2. [Conda Installation](#2-conda-installation)
   3. [Docker Installation](#3-docker-installation)
- [Usage Guide](#usage-guide)
   1. [Download](#1-download)
   2. [Build](#2-build)
   3. [Download and Build](#3-download-and-build)
   4. [Classify](#4-classify)
   5. [Profile](#5-profile)
- [Input/Output Formats](#inputoutput-formats)
- [Performance Optimization](#performance-optimization)
- [FAQ](#faq)
- [Contact & Support](#contact--support)

---

## Project Overview

**Chimera** is a versatile **metagenomic classification tool** developed by **Qinzhong Tian**, designed to simplify and accelerate the process of analyzing large-scale metagenomic datasets. Chimera integrates efficient algorithms and user-friendly features to deliver fast, accurate, and scalable metagenomic classification. 

The current version (1.1) introduces **abundance analysis**, along with the calculation of **Shannon** and **Simpson indices** for diversity assessment. Additionally, it enhances classification accuracy by incorporating the **LCA (Lowest Common Ancestor)** algorithm. Previous innovations, such as the **interleaved cuckoo filter** for constructing classification databases, ensure rapid processing and high accuracy, maintaining Chimera’s cutting-edge performance.

### 🔍 Interactive NCBI Dataset Downloads

One of Chimera’s standout features is its **interactive data downloading** capability from NCBI databases. Users can easily download and process large metagenomic datasets **within the Chimera environment**. The tool automatically handles preprocessing of downloaded datasets, streamlining the workflow from data acquisition to database construction. 

Chimera offers flexibility by supporting **custom parameter configurations**, while also providing **default settings** for users seeking a simpler setup.

### ⚡ Fast and Accurate Species Classification

Chimera is optimized for both **speed and scalability**. The classification engine is **multi-threaded**, making it highly effective at processing large datasets in a short time. The underlying **interleaved cuckoo filter** technology from version 1.0 accelerates database construction and classification, while the newly added **LCA algorithm** in version 1.1 further improves classification accuracy by resolving ambiguous taxonomic assignments.

Supported input formats include:
- Standard formats: **FASTA**, **FASTQ**
- Compressed formats: **.gz**, **.bz2**
- **Paired-end reads** for more complex data inputs

### 📊 Abundance and Diversity Analysis

Chimera version 1.1 introduces advanced **abundance analysis**, allowing users to calculate the relative abundance of taxa across multiple levels of classification. Additionally, the tool computes the **Shannon index** and **Simpson index**, providing valuable insights into species diversity and evenness within a community.

### 📊 Integrated Krona Visualization

Chimera comes with built-in **Krona integration** for visualizing taxonomic classification results. Using the **profile function**, users can easily convert their classification data into interactive **Krona charts**, allowing for **intuitive exploration** of metagenomic data. 

### 🔄 Continual Updates and Customization

Chimera is under active development, with plans for **regular updates** to introduce new features and improvements. While Chimera offers **extensive customization** for advanced users, its **default settings** ensure a simple and accessible experience for beginners. This balance makes Chimera a tool suitable for both experienced bioinformaticians and those new to metagenomic analysis.





---

## Installation

Chimera offers three installation methods: building from source, Conda installation, and Docker. We recommend using **Conda** or **source installation** for optimal performance, as Docker might introduce some overhead and reduce speed.

### 1. Source Installation

For users who prefer to build Chimera from source, here are the detailed steps. This method requires installing necessary dependencies and building Chimera manually.

#### Prerequisites

Before building Chimera, ensure you have the following dependencies installed:

- **Ubuntu 22.04** (or equivalent Linux distribution)
- **Python 3.8** (required, installable via PPA for older distributions)
- **CMake** (for compiling C++ components)
- **Krona Tools** (for visualization)
- **Essential build tools** (e.g., GCC, Make)

#### Steps

1. **Install dependencies**:

   First, update the package list and install the required build tools, Python 3.8, and other necessary libraries:

   ```bash
   sudo apt-get update
   sudo apt-get install -y software-properties-common
   sudo add-apt-repository ppa:deadsnakes/ppa
   sudo apt-get update
   sudo apt-get install -y python3.8 python3.8-dev python3.8-distutils python3-pip build-essential cmake git libbz2-dev zlib1g-dev libgcc-11-dev libstdc++-11-dev openssl libssl-dev wget bc parallel locales
   ```

2. **Install Python libraries**:

   Upgrade `pip` and install the required Python packages:

   ```bash
   python3.8 -m pip install --upgrade pip
   python3.8 -m pip install pandas multitax
   ```

3. **Install Krona Tools**:

   Download and install **Krona Tools** for visualizing classification results:

   ```bash
   wget https://github.com/marbl/Krona/releases/download/v2.8.1/KronaTools-2.8.1.tar -O /tmp/KronaTools.tar
   tar -xvf /tmp/KronaTools.tar -C /opt/
   sudo mkdir -p /opt/krona
   sudo chmod +x /opt/KronaTools-2.8.1/install.pl
   sudo /opt/KronaTools-2.8.1/install.pl --prefix /opt/krona
   sudo ln -s /opt/krona/bin/* /usr/local/bin/
   sudo ln -sf /opt/KronaTools-2.8.1/scripts/ImportText.pl /opt/krona/bin/ktImportText
   ```

4. **Clone the Chimera repository**:

   Clone the Chimera source code repository:

   ```bash
   git clone https://github.com/MalabZ/Chimera.git
   cd Chimera
   ```

5. **Build the project**:

   Create a build directory, compile Chimera, and install it:

   ```bash
   mkdir build
   cd build
   cmake ..
   make
   sudo make install
   cd ..
   ```

6. **Run Chimera**:

   After installation, you can run Chimera from the source directory using Python:

   ```bash
   python3.8 -m chimera -v
   ```

Alternatively, after installing, you can run Chimera globally using the installed `chimera.py`:

   ```bash
   chimera.py -v
   ```

This completes the source installation. You should now be able to use Chimera for metagenomic classification tasks.

### 2. Conda Installation

To install Chimera via Conda, follow these steps:

1. Create a new Conda environment with Python 3.8:
   ```bash
   conda create -n chimera python=3.8
   conda activate chimera
   ```

2. Install Chimera from the **malab** channel:
   ```bash
   conda install chimera -c malab
   ```

This method automatically resolves dependencies and is the simplest way to get Chimera running.

### 3. Docker Installation

For users preferring Docker, you can use the following commands to install and run Chimera. Note that Docker might introduce performance overhead, so for best speed, consider Conda or source installation.

1. Pull the Docker image:
   ```bash
   docker pull tianqinzhong/chimera
   ```

2. Run a test to check the installation:
   ```bash
   docker run -it --rm -v "$(pwd):/app/data" tianqinzhong/chimera -v
   ```

3. Run Chimera in Docker:
   ```bash
   docker run -it --rm -v "$(pwd):/app/data" tianqinzhong/chimera command
   ```

Replace `command` with the specific Chimera command you want to execute, such as running an analysis or building a database.


---

## Usage Guide

Chimera provides five main functions to facilitate metagenomic data processing: `download`, `build`, `download_and_build`, `classify`, and `profile`. Below are brief descriptions of each function and example usage.

### 1. Download

The `download` function allows users to fetch datasets from NCBI or other sources. Running `chimera download` enters an **interactive mode** where users can specify datasets to download.

**Example:**
```bash
chimera download
```
This command starts the interactive session for dataset downloading.

### 2. Build

The `build` function is used to construct a classification database from the downloaded datasets. It requires specifying the **input file** (usually `target.tsv` located in the downloaded folder) and allows customization of other parameters, though most parameters have sensible defaults.

**Available Parameters:**
- `-i` or `--input` (required): Input file (e.g., `target.tsv`).
- `-o` or `--output`: Output database file name (default: `ChimeraDB`).
- `-m` or `--mode`: Building mode (default: `default`).
- `-k` or `--kmer`: K-mer size (default: `19`).
- `-w` or `--window`: Window size (default: `31`).
- `-l` or `--min-length`: Minimum sequence length (default: `0`).
- `-t` or `--threads`: Number of threads (default: `32`).
- `--load-factor`: Loading ratio of ICF (default: `0.95`).
- `-q` or `--quiet`: Quiet mode.

**Example:**
```bash
chimera build -i data/target.tsv -o ChimeraDB
```
This command builds a classification database from the `target.tsv` file located in the `data/` directory and outputs the database as a single file named `ChimeraDB`.

### 3. Download and Build

The `download_and_build` function simplifies the process by downloading the required dataset and immediately building a classification database. Unlike `build`, this function does not require an input file, and the default output database name is `ChimeraDB`, though it can be customized.

**Available Parameters:**
- `-o` or `--output`: Output database file name (default: `ChimeraDB`).
- Other parameters are the same as those in the `build` function.

**Example:**
```bash
chimera download_and_build -o ChimeraDB
```
This command downloads the necessary data and directly constructs the classification database, outputting it as a single file named `ChimeraDB`.

### 4. Classify

The `classify` function allows users to perform taxonomic classification on single or paired input sequence files using the specified classification database. It supports multiple files for both single-end and paired-end reads. For paired-end reads, the number of input files must be even.

**Available Parameters:**
- `-i` or `--single`: Input files for classification (supports multiple files).
- `-p` or `--paired`: Paired input files for classification (supports multiple paired files, must be an even number).
- `-o` or `--output`: Output file name for classification results (default: `ChimeraClassify`).
- `-d` or `--database` (required): The classification database file (e.g., `ChimeraDB`).
- `-s` or `--shot-threshold`: Shot threshold for classification accuracy (default: `0.7`).
- `-t` or `--threads`: Number of threads to use during classification (default: `32`).
- `-m` or `--mode`: Classification mode, either `fast` or `normal` (default: `fast`).
- `-l` or `--lca`: Use the LCA algorithm for classification.
- `-T` or `--tax-file`: Taxonomy file for LCA classification (required if `--lca` is used, default: `tax.info` from the downloaded dataset).
- `-b` or `--batch-size`: Batch size for processing sequences (default: `400`).
- `-q` or `--quiet`: Suppresses verbose output if specified.

**Examples:**
For single-end input files:
```bash
chimera classify -i input1.fasta input2.fasta -o results.txt -d ChimeraDB
```
This command classifies the sequences in `input1.fasta` and `input2.fasta` using the `ChimeraDB` database, and saves the results to `results.txt`.

For paired-end input files (must be an even number):
```bash
chimera classify -p paired1_1.fasta paired1_2.fasta paired2_1.fasta paired2_2.fasta -o results.txt -d ChimeraDB
```
This command classifies the paired-end sequences using the `ChimeraDB` database, saving the output to `results.txt`.

For LCA-based classification with a custom taxonomy file:
```bash
chimera classify -i input.fasta -d ChimeraDB -l --tax-file tax.info -o results.txt
```
This command uses the LCA algorithm with a specified taxonomy file (`tax.info`) to classify the sequences in `input.fasta` and outputs the results to `results.txt`.

### 5. Profile

The `profile` function generates a taxonomic profile from the classification results. By default, it calculates the abundance, Shannon index, and Simpson index at different taxonomic levels (e.g., kingdom, phylum, class, order, family, genus, and species). Additionally, the `-k` option can be used to generate a **Krona chart** for interactive visualization.

**Available Parameters:**
- `-i` or `--input` (required): Input file(s) containing classification results.
- `-o` or `--output`: Output file name for the profile (default: `ChimeraProfile`).
- `-k` or `--krona`: Generate a Krona chart for interactive visualization.

By default, Chimera calculates:
- **Taxonomic abundance** at different levels.
- **Shannon index**: A measure of diversity within a community.
- **Simpson index**: A measure of dominance in the community.

**Example:**
```bash
chimera profile -i ChimeraClassify.tsv
```
This command generates a taxonomic profile with abundance, Shannon index, and Simpson index from the classification results in `ChimeraClassify.tsv`.

To generate a Krona chart:
```bash
chimera profile -i ChimeraClassify.tsv -k
```
This command generates both the taxonomic profile and a Krona chart (`ChimeraProfile.html`) for visualizing the results.

---


## Input/Output Formats

Chimera supports a variety of input and output formats to handle sequence data and results from classification tasks. This section provides details on the accepted formats for both input files and output results.

### Database Input and Construction

To construct a classification database, Chimera requires a specific input format for the database construction process.

- **Database input (`target.tsv`)**: To build a classification database, Chimera requires a `target.tsv` file. Each line in this file should contain the path to a species-specific FASTA file and the corresponding taxonomic ID (taxid), separated by a tab (`\t`).

**Example of `target.tsv` format**:
```
/path/to/species1.fasta   12345
/path/to/species2.fasta   67890
```

Once the database is constructed, Chimera stores it as a **binary file** using the **Cereal library**. This binary format allows for efficient storage and quick access during classification tasks.

**Example:**
```bash
chimera build -i target.tsv -o ChimeraDB
```
This command constructs a database from the `target.tsv` file and outputs a binary file `ChimeraDB` using the Cereal library.

### Classification Input and Output

Chimera accepts various sequence file formats for classification and generates results in a tab-separated values (TSV) format.

- **Input formats**:
  - **FASTA/FASTQ**: Chimera supports plain text FASTA and FASTQ files, as well as compressed `.gz` and `.bz2` formats.
  - **Paired-end reads**: When providing paired-end reads, an even number of input files must be specified with the `--paired` option.

**Example of valid input files**:
- `input.fasta`
- `input.fasta.gz`
- `input1.fastq`, `input2.fastq` (paired)

**Classification output format** (`classify` function):
- **TSV (Tab-Separated Values)**: Classification results are written to a TSV file. Each line includes the sequence identifier, the taxid with the highest hit count, followed by other taxids that meet the threshold.

    **Format**:
    ```
    sequence_identifier   highest_hit_taxid:hit_count   ...other_taxid:hit_count_above_threshold
    ```
    In **fast mode**, only the highest hit taxid and hit count are reported. Example:
    ```
    seq1    12345:10   67890:5
    seq2    12345:8
    ```
    If **LCA mode** is selected, the taxid classified using LCA will be represented as `taxid:0`. This indicates that the LCA algorithm was applied for classification.

**Example:**
```bash
chimera classify -i input.fasta -d ChimeraDB -o results.tsv
```
This command classifies the sequences in `input.fasta` using the `ChimeraDB` database and outputs the results in `results.tsv`.

### Profiling Output

Chimera’s `profile` function generates a detailed taxonomic profile at various levels (e.g., superkingdom, clade, phylum, class, order, family, genus, species). The output includes the count, relative abundance, Shannon index, and Simpson index for each taxonomic level. The output is presented in a tabular format, making it easy to interpret and analyze.

**Output format:**
The output is divided by taxonomic levels, and each section contains the following columns:
- **Level**: The taxonomic level (e.g., superkingdom, clade, phylum).
- **Taxon**: The name of the taxon at the specified level.
- **Count**: The number of sequences classified under that taxon.
- **Relative Abundance (%)**: The percentage of sequences relative to the total.
- **Shannon Index**: A measure of diversity.
- **Simpson Index**: A measure of dominance.

**Example Output:**

```
Level   Taxon             Count   Relative Abundance (%)  Shannon Index  Simpson Index

## Superkingdom Level ##
superkingdom    Archaea    110671   99.43   0.0353   0.0114
superkingdom    unclassified   639   0.57   0.0353   0.0114

## Clade Level ##
clade   TACK group    110671   99.43   0.0353   0.0114
clade   unclassified   639   0.57   0.0353   0.0114

## Phylum Level ##
phylum   Thermoproteota   110671   99.43   0.0353   0.0114
phylum   unclassified   639   0.57   0.0353   0.0114

## Class Level ##
class   Thermoprotei   110671   99.43   0.0353   0.0114
class   unclassified   639   0.57   0.0353   0.0114

## Order Level ##
order   Desulfurococcales   110671   99.43   0.0353   0.0114
order   unclassified   639   0.57   0.0353   0.0114

## Family Level ##
family   Desulfurococcaceae   110671   99.43   0.0353   0.0114
family   unclassified   639   0.57   0.0353   0.0114

## Genus Level ##
genus   Aeropyrum   110669   99.42   0.0356   0.0115
genus   unclassified   639   0.57   0.0356   0.0115
genus   Staphylothermus   2   0.00   0.0356   0.0115

## Species Level ##
species   Aeropyrum pernix   110610   99.37   0.0401   0.0125
species   unclassified   639   0.57   0.0401   0.0125
species   Aeropyrum camini   59   0.05   0.0401   0.0125
species   Staphylothermus hellenicus   2   0.00   0.0401   0.0125
```

This output provides detailed information on the distribution and diversity of sequences across various taxonomic levels.

**Krona Chart Option**:
Additionally, you can generate a Krona chart for interactive visualization using the `-k` option.

**Example:**
```bash
chimera profile -i results.tsv -o krona_chart -k
```
This command generates a `krona_chart.html` file from the classification results, which can be opened for visualizing the taxonomic profile.

You can see an example of Krona chart visualization here: [Krona example chart](https://telatin.github.io/microbiome-bioinformatics/data/krona/krona-test.html).

---


## Performance Optimization

Chimera is designed to handle large-scale metagenomic data efficiently, but performance can vary based on system configuration and dataset size. Below are some tips and recommendations for optimizing Chimera's performance during database construction, classification, and profiling.

### 1. Utilize Multi-threading

Chimera supports multi-threading, which can significantly reduce the time required for database construction and classification. You can control the number of threads using the `-t` or `--threads` parameter. The default is set to `32`, but you can adjust it based on the available CPU cores on your system.

**Example:**
```bash
chimera classify -i input.fasta -d ChimeraDB -o results.tsv -t 64
```
In this example, `64` threads are used, which can dramatically improve speed if your system has sufficient cores.

**Tip**: Set the number of threads to match or slightly exceed the number of physical CPU cores for optimal performance.

### 2. Adjust Batch Size

For classification tasks, Chimera processes sequences in batches. By default, the batch size is set to `400` sequences. Increasing the batch size can improve performance, especially when processing large datasets, as it reduces the overhead of repeatedly loading data.

You can adjust the batch size using the `-b` or `--batch-size` parameter.

**Example:**
```bash
chimera classify -i input.fasta -d ChimeraDB -o results.tsv -b 1000
```
Increasing the batch size to `1000` sequences can result in faster processing times, but be mindful of system memory limitations when dealing with very large batch sizes.

### 3. Optimize Database Construction Parameters

When building a classification database, you can adjust several parameters to optimize the process for your dataset:

- **K-mer size (`-k`)**: The default k-mer size is `19`, but you can adjust it based on the nature of your data. Smaller k-mer sizes might increase sensitivity but could also introduce more noise.
- **Window size (`-w`)**: Increasing the window size can reduce false positives but might slow down the construction process. The default is `31`, which works well for most datasets.
- **Minimum sequence length (`-l`)**: If your dataset contains very short sequences, consider adjusting the minimum sequence length to exclude them from the analysis. This can save processing time and improve accuracy.

**Example:**
```bash
chimera build -i target.tsv -o ChimeraDB -k 21 -w 35 -l 100
```
In this example, the k-mer size is increased to `21`, the window size to `35`, and sequences shorter than `100` base pairs are excluded.

### 4. Use Appropriate Classification Mode

Chimera offers two classification modes: `fast` and `normal`. The `fast` mode prioritizes speed by only reporting the top hit, whereas `normal` mode provides a more comprehensive result by including all taxids above the threshold.

- **Fast Mode**: Use `fast` mode (`-m fast`) when speed is a priority, and you only need the top classification hit.
- **Normal Mode**: Use `normal` mode (`-m normal`) for more detailed results, but be prepared for longer processing times.

**Example:**
```bash
chimera classify -i input.fasta -d ChimeraDB -o results.tsv -m fast
```
This example runs the classification in `fast` mode, optimizing for speed.

### 5. Manage Memory and Disk I/O

Large datasets can be memory-intensive, especially during database construction and classification. To avoid memory-related issues or performance bottlenecks:

- **Ensure sufficient RAM**: For large datasets, having more RAM allows Chimera to load and process data more efficiently.
- **Use SSDs**: If possible, store input files and the Chimera database on SSDs rather than HDDs. This can significantly reduce disk I/O bottlenecks and improve overall performance.

### 6. Adjust the Load Factor

During database construction, the **load factor** (`--load-factor`) controls the fill ratio of the internal interleaved cuckoo filter. The default value is `0.95`, but you can lower it to improve query performance at the expense of slightly larger database size.

**Example:**
```bash
chimera build -i target.tsv -o ChimeraDB --load-factor 0.85
```
In this example, a load factor of `0.85` is used, which can speed up classification queries at the cost of increasing the database size.

---

## FAQ

### 1. What is the difference between "fast" and "normal" classification modes?

- **Fast Mode**: In `fast` mode (`-m fast`), Chimera only reports the taxid with the highest hit count, making the classification process quicker but less detailed.
- **Normal Mode**: In `normal` mode (`-m normal`), Chimera reports all taxids that meet the threshold. This mode provides more comprehensive results but is slower.

### 2. What should I do if the database download is interrupted or fails?

If your database download is interrupted or fails during the interactive mode, you can choose to re-download the incomplete or failed data by enabling the fix-only mode. When prompted, type `y` to proceed.

**Example prompt:**
```
Enable fix-only mode (re-download incomplete or failed data) [y/N]:
```
Enter `y` to resume and fix the download.

### 3. What are the system requirements for running Chimera?

- **CPU**: Multi-core processors are recommended. Chimera supports multi-threading.
- **RAM**: At least 16GB of RAM is recommended, and larger datasets may require 64GB or more.
- **Storage**: SSDs are recommended for better performance.
- **Operating System**: Chimera is primarily tested on Linux (e.g., Ubuntu 20.04 or 22.04) but also supports other systems via Docker.

### 4. Can I use Chimera without Python?

Yes, you can use Chimera without relying on Python. Python is primarily used to provide functionality for downloading datasets and generating profiles. If you have built Chimera using Conda, you can simply use the following command to view available options:

```bash
Chimera -h
```

For source code builds, the `Chimera` executable is generated directly and can be run without Python:

```bash
./Chimera -h
```

For Docker, Chimera is the default entry point. To skip the default `chimera` command and access the Docker container's shell, use the following command:

```bash
docker run -it --rm -v "$(pwd):/app/data" --entrypoint Chimera tianqinzhong/chimera -h
```

### 5. Where is the taxfile required for LCA, and how do I interpret LCA results?

If you are using Chimera's built-in `download` function, the `taxfile` is located in the downloaded dataset folder as `tax.info`. For custom datasets, you will need to manually create a `taxfile` in a specific format. Each line of the file represents a taxonomic rank and includes the following fields, separated by tabs:

```
<taxid>   <parent taxid>   <rank>   <name>
```

- **taxid**: The unique identifier for the taxonomic entity.
- **parent taxid**: The taxid of the parent taxon in the hierarchy.
- **rank**: The taxonomic rank (e.g., species, genus, family, etc.).
- **name**: The scientific name of the taxon.

For example:
```
1       0           no rank        root
2157    131567      superkingdom   Archaea
2158    183925      order          Methanobacteriales
2159    2158        family         Methanobacteriaceae
2160    2159        genus          Methanobacterium
2162    2160        species        Methanobacterium formicicum
```

This example defines a taxonomic hierarchy starting from `root` (no rank) down to the species **Methanobacterium formicicum**.

- **taxid 1** is the root of the hierarchy with no parent (`parent taxid = 0`).
- **taxid 2157** represents the **Archaea** superkingdom, which belongs to the parent taxon **131567**.
- Similarly, **taxid 2162** represents the species **Methanobacterium formicicum**, which is a descendant of the genus **Methanobacterium** (`taxid 2160`).

### Interpreting LCA Results

In the classification output, any result classified using the LCA algorithm will be shown as `taxid:0`. This indicates that the Lowest Common Ancestor (LCA) method was applied, and the classification could not be resolved to a more specific taxonomic level.



If you have further questions, feel free to ask by opening an issue on our [GitHub repository](https://github.com/LoadStar822/Chimera/issues).



---

## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

---

## Contact & Support
For any questions or support, feel free to reach out to us:
- **Website**: [MalabZ](http://lab.malab.cn/~cjt/MSA/)
- **Email**: tianqinzhong@qq.com
