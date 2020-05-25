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
# needs to be run only once
~$ conda create --name env_name --file requirements.txt

# activate environment
~$ conda activate env_name

... run analysis ...
# deactivate environment when done
~$ conda deactivate
```

### Usage

`bard` does not work directly on the raw read sequence files. They need to be preprocessed (removing rRNAs) and mapped to a reference genome first, using either `bowtie2/tophat` or some similar software of your choice. You'll also need the corresponding annotation file and a multiple fasta file containing all the CDS sequences. Please be aware that the gene names in the GTF/GFF (annotation) file need to exactly match those in the CDS file. If you are working on yeast, the chromome IDs in your reference genome should match with those in the GTF.

It is preferable to download these organism specific files from the same source, such as `Ensembl`.

`TL; DR`: You need to have the following files available for input:
* An annotation file (either GTF/GFFv3)
* An alignment file in BAM format (remove unmapped reads and index it first, using `samtools` or alike)
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
The following is a short description of how to run the analysis for *S. cerevisiae*. The read sequences were obtained from the library `SRR1042864`. The raw sequences were preprocessed and mapped to the *S.cerevisiae* `S288C_R64-1-1_20110203` reference genome according to the protocol detailed in `10.1038/nprot.2012.086`. All organism specific files were obtained from `yeastgenome.org`.

The configuration file we'd use in this case:

```json
{
    "coding_sequence_path": "/full/path/to/cds.fa",
    "coding_sequence_format": "fasta",
    "annotation_file_path": "/full/path/to/annotation.gff",
    "annotation_feature_tag": "protein_id",
    "bam_file_path": "/full/path/to/alignment.bam",
    "read_offset_terminal": "three_prime",
    "coverage_cutoff": 40,
    "coverage_metric": "reads_per_nt",
    "will_ignore_overlaps": true,
    "peak_scan_range": [5, 30],
    "use_readlengths": [26, 27, 28, 29],
    "gene_list_file": "",
    "gene_list_action": "",
    "genes_overlap_exception": ""
}
```
Briefly, the options are:
* `annotation_feature_tag`:  The `protein_id` tag in the **9th** column of the GTF file uniquely identifies the gene (feature) names
* `read_offset_terminal`: Determine the E/P/A site offsets from the `3'` ends of the reads
* `coverage_cutoff`: Ignore all genes with a read coverage below this numerical threshold (40)
* `coverage_metric`: Calculate the read coverage as average number of reads mapping per nucleotide
* `will_ignore_overlaps`: Ignore all overlapping genes
* `peak_scan_range`: Check for the metagene initiation peak within this range of nucleotides w.r.t the annotated start site.
* `use_readlengths`: Ignore all reads which are not of the given sizes (in nt).


The results will be written to a separate folder in your working directory. It has the following structure:
```shell
# For illustrative purposes
# File names will differ
.
├── logfile.txt
├── configuration_used.json
│
├── Data
│   ├── data_file_1.json
│   ├── data_file_2.json
│   └── data_file_3.json
│            ... and so on
└── Plots
    ├── plot_1.svg
    ├── plot_2.svg
    └── plot_3.svg
             ... and so on
```

The data files are JSON files which you can load into a Python or R environment. For the above analysis, we would get the folllowing results:

##### Initiation peak
<img src=/images/initiation.png width=400 height=150> <img src=/images/start_framing.png width=400 height=150>

##### Reading frame
<img src=/images/reading_frame.png width=500 height=300> <img src=/images/meta_rpf_framing.png width=350 height=400>

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




## Contact

Queries can be directed to

## License
GPLv3
