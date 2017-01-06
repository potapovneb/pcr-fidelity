# Examining Sources of Error in PCR by Single-Molecule Sequencing

**Vladimir Potapov and Jennifer L. Ong**

## Description
This repository provides a set of scripts for a PacBio single-molecule sequencing assay used to comprehensively catalog the different types of errors introduced during PCR. The full details of the method are available in the following publication:

Potapov V. & Ong JL. Examining Sources of Error in PCR by Single-Molecule Sequencing. PLOS ONE. 2016. doi:10.1371/journal.pone.0169774. ([View article](http://dx.doi.org/10.1371/journal.pone.0169774))

## Usage
The analysis workflow consists of several consecutive steps as outlined below. The correspondng scripts can be found in `scripts/` directory.

1. Generating strand-specific consensus sequences (`ccs2-consensus.sh`).

2. Mapping consensus sequences to a reference sequence and extracting mutations (`ccs2-mutation.sh`).

3. Filtering and summarizing mutations (`ccs2-chimeric.sh`,`ccs2-summary.sh`)

4. Additional scripts to analyze template switching and PCR-mediated recombination can be found in `extra/` directory.

The included [workflow](workflow.md) provides detailed instructions for data analysis conducted in the original publication [(1)](#ref1).

## Requirements <a name="requirements"></a>
* [SMRT Link](https://github.com/PacificBiosciences/SMRT-Link). Scripts rely on `bax2bam`, `ccs` and `blasr` command-line utilities, which are developed by Pacific Biosciences, Inc. and distributed as part of the SMRT Link software.
* [SAMtools](http://samtools.sourceforge.net/). Manipulating BAM files. 
* [BWA](http://bio-bwa.sourceforge.net/). Detecting chimeric reads.
* [P7ZIP](http://p7zip.sourceforge.net/). Compressing output files to minimize disk usage.
* [The R Project for Statistical Computing](https://www.r-project.org/). Extracting and tabulating data.

## Citations<a name="ref1"></a>
1. Potapov V. & Ong JL. Examining Sources of Error in PCR by Single-Molecule Sequencing. PLOS ONE. 2016. doi:10.1371/journal.pone.0169774. ([View article](http://dx.doi.org/10.1371/journal.pone.0169774))
