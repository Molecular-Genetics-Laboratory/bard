#  bard
**bard** is a command-line tool for analyzing ribosome profiling datasets. It can calculate E/P/A site offsets, metagene profiles, detect pauses, and more. The current version works for prokaryotes and yeast.


## Requirements
* Python 3.x (>= 3.7)
* Biopython (>= 1.76)
* Numpy (>= 1.18.1)
* Seaborn (>= 0.10.0)
* Matplotlib (>= 3.1.3)
* Pysam (>= 0.15.3)
* Linux / macOS
* At least 2GB of (usable) RAM


## Installation and usage
After ensuring that all dependencies are satisfied, you can directly download and run the `bard.py` script.
### Installing dependencies
If you don't have the required libraries already installed, you can do so via `pip` using the `requirements.txt` file:

```shell
~$ pip install requirements.txt
```

..or you could set up a virtual environment using `conda`:

```shell
# create virtual environment
~$ conda create --name env_name --file requirements.txt

# activate virtual environment
~$ conda activate env_name

... run analysis ...
# deactivate environment when done
~$ conda deactivate
```

### Usage

`bard` does not work directly on the raw read sequence files. They need to be preprocessed (removing rRNAs) and mapped to a reference genome first, using either `bowtie2/tophat` or some similar software of your choice. You'll also need the corresponding annotation file and a multiple fasta file containing all the CDS sequences. Please be aware that the gene names in the GTF/GFF (annotation) file need to exactly match those in the CDS file. Also, if you are working on yeast, the chromome IDs in your reference genome should match with those in the GTF.

It is preferable to download these organism specific files from the same source, such as `Ensembl`.

`TL; DR`: You need to have the following files available for input:
* An annotation file (either GTF/GFFv3)
* An alignment file in BAM format (please remove unmapped reads and index it first, using `samtools` or alike)
* A multi-fasta file containing the CDS sequences for each gene

#### First run
`bard` does not have any command-line options; using only a single configuration file instead. To create a configuration file, run the script without any arguments:

```shell
~$ python bard.py
```

This will create two files in your working directory: `bard_config_template.json`, which is an empty configuration file you can edit, and a short help file `bard_config_help.txt`, which describes the available options. Please go through the help file for full details.

The config file looks something like this:

```json
{
    "coding_sequence_path": "/home/user/path/to/cds.fa",
    "coding_sequence_format": "fasta",
    "annotation_file_path": "/home/user/path/to/annotation.gff",
    "annotation_feature_tag": "ID",
    "bam_file_path": "/home/user/path/to/alignment.bam",
    "read_offset_terminal": "five_prime",
    "coverage_cutoff": 40,
    "coverage_metric": "reads_per_nt",
    "will_ignore_overlaps": true,
    "peak_scan_range": [-25, -5],
    "use_readlengths": [26,27,28,29,30,31,32,33],
    "gene_list_file": "/home/user/path/to/gene_list.txt",
    "gene_list_action": "include_only",
    "genes_overlap_exception": "/home/user/path/to/overlapping_genes_list.txt"
}
```
Once you have edited the configuration file to your requirements, your can run `bard` like so:

```shell
# time for a short coffee break :)
~$ python bard.py config.json
```

#### Example (for *S. cerevisiae*)
A short tutorial on how to do this for yeast.

```shell
.
├── logfile.txt
├── configuration_used.json
├── Data
│   ├── Annotation.json
│   ├── Asite_Offset.json
│   ├── Asymmetry_Scores.json
│   ├── Reading_Frame_Fractions.json
│   ├── Data_for_Termination_Peak.json
│   ├── Esite_Offset.json
│   ├── Gene_Coverages.json
│   ├── Reading_Frame_Rawvalues.json
│   ├── High_Coverage_Genes.json
│   ├── Read_Coverage_Profiles.json
│   ├── Metagenes_per_read_length_without_offset.json
│   ├── Metagenes_per_read_length_with_offset.json
│   ├── Amino_Acid_Pause_Scores.json
│   ├── Codon_Pause_Scores.json.json
│   ├── Density_vector_per_Gene.json
│   ├── Metagene_of_positional_pause_scores.json
│   ├── Positional_Pause_scores_per_Gene.json
│   ├── Psite_offset.json
│   ├── Read_lengths.json
│   ├── Data_for_start_peak.json
│   └── Terminal_Base_Fractions.json
└── Plots
    ├── Base_Fraction_ThreePrime_End.tif
    ├── Base_Fraction_FivePrime_End.tif
    ├── Asymmetry_Scores.tif
    ├── Framing_Metagene_per_read_length_no_offset.tif
    ├── Framing_Metagene_per_read_length_offset.tiff
    ├── Metagene_no_offset.tiff
    ├── Metagene_with_offset.tiff
    ├── Pause_score_per_aminoacid.tiff
    ├── Pause_score_per_codon.tiff
    ├── Positional_pause.tiff
    ├── Global_reading_frame.tiff
    ├── Read_length_histogram.tiff
    ├── Initiation_peak.tiff
    ├── Framing_initiation_peak.tiff
    └── Termination_peak.tiff

```

### Processing multiple BAM files in parallel
At present bard is single threaded. If you wish to process multiple files simultaneously, then place all the config files in a single folder, then run:

```shell
~$ for config in *.json; do echo python bard.py $config; done > jobs.txt

```

This will create a jobs.txt file with the following contents:

```text
python bard.py config_1.py
python bard.py config_2.py

...

python bard.py config_N.py
```

Then you can use GNU Parallel to run them simultaneously:

```shell
~$ parallel < jobs.txt
```

## Configurations
Instead of having a plethora of command line switches, bard uses a JSON file for runtime configurations.

If you don't have a configuration file, do this:
```shell
~$ python bard.py
```

bard will helpfully print out the configuration template and associated help text. Let's explore the options in more detail:

* **read_length_histogram**:   true/false. Do you want to plot the historam?  

## TODO
* Add support for eukaryotic organisms
* Parallel processing of multiple BAM files without help from external programs
* Documentation


## Contact

Queries can be directed to

## License
GPLv3
