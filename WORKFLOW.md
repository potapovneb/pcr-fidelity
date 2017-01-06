# Examining Sources of Error in PCR by Single-Molecule Sequencing

**Vladimir Potapov and Jennifer L. Ong**

## General instructions

1. Obtain scripts:
   ```
   git clone https://github.com/potapovneb/pcr-fidelity /path/to/downloaded-code
   ```

2. Set `$PROJDIR` in `workflow.sh` to match the location of the top directory of the downloaded code:<a name="projdir"></a>
   ```
   export PROJDIR=/path/to/downloaded-code
   ```

3. Download raw sequencing data to `pacbio_data/` directory (SRA accession number SRP095133). Consult with `input/samples.md5` for correct location and data integrity of sequencing files.

4. Make sure that all [required tools](README.md#requirements) are available in your system. Analysis scripts use a set of
   ```
   module load samtools
   module load p7zip
   etc.
   ```
   commands to manipulate `$PATH` variable. This can be easily replaced by directly setting `$PATH` environment variable:
   ```
   export PATH=/path/to/external-tools:$PATH
   ```

5. Analysis scripts were designed to analyse a single SMRTCell at a time using Open Grid Scheduler/Grid Engine. This allows processing several SMRTCells in parallel. If OGS/GE is not available, invoking scripts through `qsub` has to be changed.

## Identifying polymerase base substitution errors in PCR

1. Build strand specific-consensus reads:
   ```shell
   ./workflow.sh input/samples.csv consensus
   ```

1. Extract mutations:
   ```shell
   ./workflow.sh input/samples.csv mutation
   ```

1. Detect chimeric reads:
   ```shell
   ./workflow.sh input/samples.csv chimeric
   ```

1. Filter and summarize extracted data:
   ```shell
   ./workflow.sh input/samples.csv summary
   ```

1. Build tables:
	```shell
	./workflow.sh input/samples.csv tabulate > results/results.csv
	cd results
	Rscript ../scripts/build_tables.R results.csv
	```

## Template switching at cruciform DNA structures in *lacZ*

1. Switch to subproject directory:
	```
	cd extra/inversions
	```

2. Set `$PROJDIR` in `workflow.sh` as [described](#projdir).

3. Run the analysis:
	```
	./workflow.sh execute
	```

4. Extract data and build tables:
	```
	./workflow.sh extract > results.csv
	```

## PCR-mediated recombination

1. Switch to subproject directory:
	```
	cd extra/pcr_mediated_recombination
	```

2. Set `$PROJDIR` in `workflow.sh` as [described](#projdir).

3. Run the analysis and build tables:
	```
	./workflow.sh > results.csv
	```
 
## Citations<a name="ref1"></a>
1. Potapov V. & Ong JL. Examining Sources of Error in PCR by Single-Molecule Sequencing. PLOS ONE. 2016. doi:10.1371/journal.pone.0169774. ([View article](http://dx.doi.org/10.1371/journal.pone.0169774))
