<p align="center">
	<img src="https://github.com/user-attachments/assets/119df30c-61ab-43a5-9b66-62df67708578" alt="KATSS Logo of a blue cat with a RNA-shaped tail." width="300">
</p>

# KATSS - K-mer Analysis Tools for Sequence and Structure

A C package containing programs to analyze RNA-protein interaction.

## Introduction

KATSS is a collection of C programs designed for the analysis of RNA sequences, particularly focusing on RNA-protein interactions. The package includes two main programs: `kstruct` and `ikke`, each tailored for specific analyses.

The tools allow you to:
* Calculate the motif of RNA-binding proteins
* Determine the secondary motifs
* Statistically determine enrichments without a control dataset
* Predict the binding preference of a protein

## Table of Contents
1. [Installation](#installation)
2. [Executable Program](#executable-programs)
3. [Example Usage](#example-usage)
4. [License](#license)
## Installation

KATSS uses CMake to install the programs. In order to compile from source, the following dependencies are required:
* ```ViennaRNA```, for RNA structure prediction
* ```zlib``` OR ```ISA-L```, to read from compressed files
* ```CMake >= 3.9.0```, the build toolchain used by katss

### Quick Start

Once you download the project, you can simply do the following:

```bash
cd katss
cmake -B build && cmake --build build
cmake --build build --target install
```

**Note**: This will install the binaries into the default system path, which is typically `/usr/local/bin`.

### User-dir Installation

One issue you might encounter with the previous installation is that it might require root privileges. In case you do
not have root privileges in your computer, you can specify the installation location to a place you have write access
to by using the ```-DCMAKE_INSTALL_PREFIX``` flag as such:
```bash
cmake -DCMAKE_INSTALL_PREFIX=/your/preferred/path -B build && cmake --build build
cmake --build build --target install
```

Ensure that the specified path is included in your system's `PATH` environment variable to execute the binaries conveniently from any location.

## Executable Programs

The KATSS package includes the following executable programs:
| Program   | Description                                                              |
| --------- | ------------------------------------------------------------------------ |
| `ikke`    | Compute enriched motifs in RNA-sequences                                 |
| `kstruct` | Predict k-mer preference in base-pair interactions in protein-bound RNA  |

To get more information about a program, run them with the `--help` or `--detailed-help` flag as such:

```bash
ikke --detailed-help
```

## Example Usage

### `kstruct`

```bash
kstruct -i control_sequences.fastq.gz -b bound_sequences.fastq.gz -o output.csv -k 3
```

### `ikke`

```bash
ikke -t test_seqs.fastq.gz -c ctrl_seqs.fastq.gz -o output -dt --kmer=6 --iterations=10
```

## License

This project is licensed under GNU General Public License v3.0.
See the [LICENSE](LICENSE) file for full details.
