# Chimera

**Ultrafast and memory-efficient database construction for high-accuracy metagenomic taxonomic classification and profiling.**

[![License](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-linux--64-lightgrey)](https://www.linux.org/)
[![Conda](https://img.shields.io/conda/vn/malab/chimera)](https://anaconda.org/malab/chimera)
[![Preprint](https://img.shields.io/badge/preprint-bioRxiv-orange)](https://doi.org/10.1101/2025.03.26.645388)

Chimera is a reference-database metagenomic classifier designed for fast custom database construction, accurate read-level taxonomic assignment, and abundance profiling in the same run. It supports both long and short reads, builds databases from NCBI genomes or custom references, and can apply local near-neighbor resolution for difficult strain-like assignments.

<p align="center">
  <img src="assets/chimera_usage_demo.gif" alt="Chimera build, classify, and profile demo" width="900">
</p>

<p align="center">
  <em>Database build, read classification, and profile generation with Chimera.</em>
</p>

## Highlights

- **Fast custom database construction** from downloaded NCBI genomes or a user-provided `target.tsv`.
- **High-accuracy read classification** for large reference-database workflows.
- **Profile output during classification**: `chimera classify` writes both per-read assignments and `ChimeraProfile.tsv`.
- **Interactive NCBI download wizard** for users who do not want to hand-write genome download commands.
- **Automatic NCBI taxdump handling** during build when the Python wrapper can infer or download the taxonomy data.
- **Optional local read resolution (LPC)** for near-neighbor or strain-like ambiguity when the database contains LPC data.
- **Simple user workflow**: interactive genome download, automatic taxonomy handling, minimal required parameters, and classify/profile output in one command.

## Quick Links

- [Install](#install)
- [Quick Start](#quick-start)
- [Download Genomes](#download-genomes)
- [Build a Database](#build-a-database)
- [Classify and Profile](#classify-and-profile)
- [Custom `target.tsv`](#custom-targettsv)
- [Output Files](#output-files)
- [Citation](#citation)

## Install

Conda is the recommended installation route for normal use. The Chimera package is currently published on the `malab` Anaconda channel; `conda-forge` and `bioconda` are used for dependencies.

```bash
conda create -n chimera -c malab -c conda-forge -c bioconda chimera
conda activate chimera
chimera -v
```

The prebuilt Conda binary targets Linux x86-64 machines with AVX2 support. If your machine does not support AVX2, build from source with portable SIMD settings:

```bash
git clone --recursive https://github.com/LoadStar822/Chimera.git
cd Chimera

python -m pip install -e .

cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCHIMERA_ENABLE_NATIVE_AVX2=OFF
cmake --build build -j 32

export CHIMERA_BIN="$PWD/build/Chimera"
chimera -v
```

For AVX2-capable source builds, omit `-DCHIMERA_ENABLE_NATIVE_AVX2=OFF`. Link-time optimization is disabled by default for portability and can be enabled with `-DCHIMERA_ENABLE_LTO=ON`.

## Quick Start

The Python `chimera` command is the standard user-facing entry point. It wraps download, build, output directory checks, and taxonomy handling around the native Chimera engine.

```bash
# 1. Download reference genomes from NCBI.
#    With no arguments, this opens the interactive downloader.
chimera download

# 2. Build an indexed Chimera database.
chimera build \
  -i genome_output/target.tsv \
  -o ChimeraDB \
  -t 32

# 3. Classify reads and write an abundance profile.
chimera classify \
  -i reads.fastq.gz \
  -d ChimeraDB \
  -o results \
  -t 32
```

After classification, the output directory contains:

```text
results/
  ChimeraClassify.tsv
  ChimeraProfile.tsv
```

For CAMI/OPAL-compatible benchmark output, add `--profile-cami`:

```bash
chimera classify \
  -i reads.fastq.gz \
  -d ChimeraDB \
  -o results \
  -t 32 \
  --profile-cami
```

This additionally writes `ChimeraProfile.cami.tsv`.

## Typical Workflow

| Step | Command | Main result |
| --- | --- | --- |
| Download genomes | `chimera download` | `genome_output/target.tsv`, `taxdump/`, genome files |
| Build database | `chimera build -i genome_output/target.tsv -o ChimeraDB` | `ChimeraDB/` |
| Classify reads | `chimera classify -i reads.fastq.gz -d ChimeraDB -o results` | `ChimeraClassify.tsv`, `ChimeraProfile.tsv` |
| Benchmark export | add `--profile-cami` | `ChimeraProfile.cami.tsv` |

## Download Genomes

Run:

```bash
chimera download
```

Without additional arguments, Chimera starts an interactive NCBI download wizard. It asks for the source database, organism group, assembly filters, output directory, thread count, and optional taxid filters.

The downloader writes a build-ready directory:

```text
genome_output/
  target.tsv
  tax.info
  taxonomy.meta
  taxdump/
  2026-.../
    files/
      ...
```

For scripted downloads, pass the same choices as command-line options:

```bash
chimera download \
  -d refseq \
  -g archaea,bacteria \
  -l "complete genome" \
  -c "reference genome" \
  -f genomic.fna.gz \
  -o genome_output \
  -t 16
```

Useful options:

| Option | Meaning |
| --- | --- |
| `-d refseq` | Use RefSeq. |
| `-g archaea,bacteria` | Select organism groups. |
| `-T 562,623` | Restrict to NCBI taxids. |
| `-A 100` | Cap assemblies for a small test. |
| `-k` | Dry run. |

The documented workflow in this README uses NCBI taxonomy.

## Build a Database

Build from downloader output:

```bash
chimera build \
  -i genome_output/target.tsv \
  -o ChimeraDB \
  -t 32
```

`ChimeraDB` is a directory. Pass the directory itself to `classify`:

```bash
chimera classify -i reads.fastq.gz -d ChimeraDB -o results
```

When `target.tsv` comes from `chimera download`, the wrapper finds the adjacent `taxdump/`. If NCBI taxdump is missing, the wrapper can download and verify it automatically before invoking the native build.

### Local Read Resolution

For NCBI databases, Chimera builds local read resolution data by default when usable taxonomy data are available. LPC is used by `classify` only when the database contains the required data and the sample contains eligible local ambiguity.

Disable LPC during build if you want a smaller database or do not want this extra build step:

```bash
chimera build \
  -i genome_output/target.tsv \
  -o ChimeraDB \
  --no-local-resolution
```

Disable LPC during classification if you want to force the main classifier path:

```bash
chimera classify \
  -i reads.fastq.gz \
  -d ChimeraDB \
  -o results \
  --no-local-resolution
```

## Custom `target.tsv`

You can build Chimera databases from your own references. Write a two-column `target.tsv` with no header:

```text
/absolute/path/to/reference_1.fna.gz    562
/absolute/path/to/reference_2.fna.gz    562
/absolute/path/to/reference_3.fna.gz    623
```

Column 1 is a FASTA/FASTQ file path. Column 2 is the taxid assigned to all sequences in that file. Tabs are recommended; any whitespace separator is accepted by the build parser. Multiple files can use the same taxid.

For NCBI builds, use NCBI taxids. If you provide a custom `target.tsv`, either place a matching `taxdump/` beside it or pass `--taxonomy-dir`:

```bash
chimera build \
  -i target.tsv \
  -o ChimeraDB \
  --taxonomy-kind ncbi \
  --taxonomy-dir /path/to/taxdump \
  -t 32
```

If `--taxonomy-dir` is omitted and the build uses NCBI taxonomy, the Python wrapper can download NCBI taxdump automatically.

## Classify and Profile

Single-end reads:

```bash
chimera classify \
  -i reads.fastq.gz \
  -d ChimeraDB \
  -o results \
  -t 32
```

Paired-end reads:

```bash
chimera classify \
  -p sample_R1.fastq.gz sample_R2.fastq.gz \
  -d ChimeraDB \
  -o results \
  -t 32
```

With the Python wrapper, `-o` is an output directory. The directory must be empty or absent; Chimera will not write into a non-empty existing output directory.

`classify` already writes `ChimeraProfile.tsv`. A separate `profile` command is only needed for auxiliary conversion from an existing aggregate table:

```bash
chimera profile -i existing_aggregate.tsv -o ChimeraProfile
chimera profile -i existing_aggregate.tsv -o ChimeraProfile -k
```

The `-k` option requests Krona output. Krona Tools are optional and are not installed by the default Chimera Conda package; install Krona separately if you need HTML Krona plots.

## Output Files

| File | Written by default | Purpose |
| --- | --- | --- |
| `ChimeraClassify.tsv` | yes | Per-read taxonomic assignments and read-level metadata. |
| `ChimeraProfile.tsv` | yes | Default abundance profile for normal Chimera use. |
| `ChimeraProfile.cami.tsv` | only with `--profile-cami` | CAMI/OPAL-compatible exchange format for benchmark tools. |

### `ChimeraClassify.tsv`

Per-read classification output. Each row reports the read identifier, selected taxonomic assignment, and optional metadata such as posterior top-k entries or rejection reasons.

Example shape:

```text
read_id    taxid:score    POST_TOPK=...
read_id    unclassified   REJECT=...   HINT=...
```

This file is for read-level inspection and downstream filtering. It is not the abundance profile, and `ChimeraProfile.tsv` should not be reconstructed by simply counting final labels.

### `ChimeraProfile.tsv`

Taxonomic abundance profile produced by `classify`.

Main columns:

| Column | Meaning |
| --- | --- |
| `rank` | Taxonomic rank of the reported row. |
| `taxid` | NCBI taxid. |
| `name` | Scientific name. |
| `relative_abundance` | Abundance fraction. |
| `percentage` | Abundance percentage. |
| `read_support` | Count-like profile support; not a direct histogram of final labels. |
| `parent_taxid`, `parent_name` | Parent lineage context. |
| `genus_taxid`, `genus_name` | Genus-level context when available. |
| `lineage_taxids`, `lineage_names` | Semicolon-separated lineage. |
| `single_source_limited` | Whether the row is supported by only one reference/source context in the profile readout. |
| `reportability` | Why the row is reportable in the default profile. |

### `ChimeraProfile.cami.tsv`

CAMI/OPAL-compatible profile output. It is written only when `--profile-cami` is requested.

## Database Layout

A Chimera database is a directory containing the core index, taxonomy metadata, and optional LPC data:

```text
ChimeraDB/
  core.imcf
  manifest.tsv
  ...
```

Internal filenames may change between releases. For reproducibility, keep the whole database directory together rather than moving only `core.imcf`.

## Performance Notes

- Put the database and reads on SSD storage when possible.
- Use `-t` to set thread count. The wrapper default is `min(os.cpu_count(), 64)`. You can pass a higher value explicitly, but more threads are not always faster.
- Large reference databases can require substantial RAM during both build and classify.
- LPC can improve difficult local assignments but adds runtime and memory cost.
- Prebuilt binaries target AVX2-capable x86-64 machines. Use a portable source build on older CPUs.

## Benchmarking

Larger tool comparisons are maintained in a separate repository:

```text
https://github.com/LoadStar822/ChimeraBenchmark
```

The benchmark repository is optional and is not required for normal download, build, or classify workflows. If you clone Chimera with `--recursive`, Git may initialize the benchmark submodule under `third_party/`; ordinary users can ignore it.

For benchmark tools that require CAMI format, run `classify` with `--profile-cami`.

## Citation

If you use Chimera, please cite:

**Chimera: Ultrafast and Memory-efficient Database Construction for High-Accuracy Taxonomic Classification in the Age of Expanding Genomic Data**

- BioRxiv: <https://www.biorxiv.org/content/10.1101/2025.03.26.645388v1>
- DOI: <https://doi.org/10.1101/2025.03.26.645388>

GitHub also reads the repository `CITATION.cff` file.

## References and Acknowledgements

Chimera builds on several open-source projects:

- [SeqAn3](https://github.com/seqan/seqan3) provides modern sequence I/O, sequence alphabets, and BGZF-aware stream support used throughout the build and classify pipeline.
- [CLI11](https://github.com/CLIUtils/CLI11) provides the C++ command-line parser.
- [cereal](https://github.com/USCiLab/cereal) and [SDSL-lite](https://github.com/simongog/sdsl-lite), distributed with SeqAn3, support serialization and succinct data-structure components used by Chimera indexes.
- [robin-hood hashing](https://github.com/martinus/robin-hood-hashing) is used for fast hash tables in classification and local-resolution code paths.
- [xxHash](https://github.com/Cyan4973/xxHash) is used for fast non-cryptographic hashing.
- [SIMDe](https://github.com/simd-everywhere/simde) provides portable SIMD interfaces used when building Chimera without AVX2-specific code.
- [strobemers](https://github.com/ksahlin/strobemers) and ideas from [strobealign](https://github.com/ksahlin/strobealign) inform Chimera's strobemer/randstrobe feature construction.
- [genome_updater](https://github.com/pirovc/genome_updater) is bundled for NCBI genome download workflows.
- [Krona Tools](https://github.com/marbl/Krona) can be used by the auxiliary profile conversion command when Krona output is requested.
- [klib](https://github.com/attractivechaos/klib) contributes low-level C utilities used in the C/C++ codebase.
- [OpenMP](https://www.openmp.org/), [zlib](https://zlib.net/), [bzip2](https://sourceware.org/bzip2/), and [OpenSSL](https://github.com/openssl/openssl) provide parallel execution, compressed file support, and checksum/download infrastructure.

We thank the maintainers of these projects and the broader open-source bioinformatics community.

## License

Chimera is released under the MIT License. See [LICENSE](LICENSE).

## Contact

- GitHub: <https://github.com/LoadStar822/Chimera>
- Author: Qinzhong Tian
- Homepage: <https://loadstar822.github.io/>
