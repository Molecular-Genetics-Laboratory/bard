#  bard
**bard** is a tool for visualization and reporting on ribosome profiling datasets. The curent version is compatible with prokaryotes and yeast.

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
Once all dependencies are satisfied, you can download and run the `bard.py` script directly.
### Installing dependencies
You can install the required libraries via `pip` using the `requirements.txt` file:

```shell
~$ pip install requirements.txt
```

..or set up a virtual environment using `conda`:

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

`bard` does not work with the raw read sequence (FASTQ) files directly. They need to be preprocessed (removing rRNAs, adapters, etc) and mapped to a reference genome first, using either `bowtie2 / tophat` or some similar tool of your choice. You'll also need the corresponding annotation file and a multiple fasta file containing all the CDS sequences. The gene names in the GTF/GFF (annotation) file need to exactly match those in the fasta file. If you are working on yeast, the chromome IDs in your reference genome (fasta) should match with those in the GTF (first column).

It is preferable to download these organism specific files from the same source, such as `Ensembl`.

So, you need to have the following files ready:
* An annotation file (either GTF/GFFv3)
* An alignment file in BAM format (indexed and unmapped reads removed)
* A multi-fasta file containing the CDS sequences for each gene

#### First run
`bard` does not have any command-line options; using only a single configuration file instead. To create a configuration file, run the script without any arguments:

```shell
~$ python bard.py
```

This will create two files in your working directory: `bard_config_template.json`, which is an empty configuration file you can edit, and a short help file `bard_config_help.txt`, which describes the available options. Please go through the help file for a full description.

The configuration file looks something like this:

```json
{
    "coding_sequence_path": "/home/user/path/to/cds.fa",
    "coding_sequence_format": "fasta",
    "annotation_file_path": "/home/user/path/to/annotation.gff",
    "annotation_feature_tag": "ID",
    "bam_file_path": "/home/user/path/to/alignment.bam",
    "check_offset_from": "five_prime",
    "coverage_cutoff": 10,
    "coverage_metric": "rpkm",
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
# coffee?
~$ python bard.py config.json
```

### Example usage
The following is a short description of how to run the script for a sample  *S. cerevisiae* dataset. The read sequences were obtained from the library `SRR1042864`. The raw sequences were preprocessed and mapped to the *S.cerevisiae* `S288C_R64-1-1_20110203` reference genome according to the protocol detailed in `doi:10.1038/nprot.2012.086`. All organism specific files were obtained from `yeastgenome.org`.

We'd use the following configuration:

```json
{
    "coding_sequence_path": "/home/user/Desktop/Yeast_S288C_CDS.fa",
    "coding_sequence_format": "fasta",
    "annotation_file_path": "/home/user/Desktop/Yeast_S288C_annotation.gff",
    "annotation_feature_tag": "protein_id",
    "bam_file_path": "/home/user/Desktop/aligned_reads.bam",
    "check_offset_from": "three_prime",
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

Option | Description
-------|------------
`bam_file_path`| Absolute path to the read alignment (BAM) file
`coding_sequence_path`| Absolute path to the multiple fasta file containing CDS sequences for each gene
 `annotation_feature_tag`|  This tag in the `9th` GTF attribute column specifies the gene names (`protein_id`)
`check_offset_from`| Determine E/P/A-site offsets from the `3'` ends of the ribosome protected fragments (reads)
`coverage_cutoff`| Ignore all genes with a read coverage below this numerical threshold (`40`)
`coverage_metric`| Calculate the read coverage as average number of reads mapping per nucleotide (`reads_per_nt`)
`will_ignore_overlaps`| Ignore all overlapping genes (`true`)
`peak_scan_range`| Check for the translation initiation signal within `5 to 30` nt downstream of annotated start site
`use_readlengths`| Only use reads of lengths `26, 27, 28, 29` nt in the analysis
`gene_list_file`| Absolute path to text file with a list of genes (`none specified`)
`gene_list_action`| What to do with the  provided gene list (`no action specified`)
`genes_overlap_exception`| Include these overlapping genes in the analysis (`none specified`)


Run the script the same way as mentioned above. The results will be written to a separate folder in your working directory with the following structure:
```shell
# For illustrative purposes
# File names will differ
<id>_Results_<bam_file_name>_<timestamp>
│
├── <id>_logfile.txt       # full runtime logs
├── <id>_runconfig.json    # the configuration in it's entirity
│
├── Data
│   ├── <id>_annotation.json       # data available for further processing
│   ├── <id>_psite_offsets.json
│   └── <id>_gene_coverages.json
│
│            ... and so on
└── Plots                         # analysis plots
    ├── <id>_pause_score_codon.svg
    ├── <id>_reading_frame.svg
    └── <id>_asymmetry_scores.svg

             ... and so on
```
`<id>` will be unique each time you run `bard`. This makes it possible to run multiple instances in parallel (see below) without the risk of overwriting existing files.

The JSON files under `Data/` can be loaded into a Python or R (using `jsonlite`, for instance) environment for further analysis.

### Resulting output
A summary of the translation process being investigated gets written out as an [HTML report.](../blob/master/examples/report.html). Some relevant plots are described below:

#### Basic diagnostics
Plots showing the overall quality of the riboseq dataset. Here, the top left plot shows the global reading frame, where ~80% of the reads (p-sites) map to `frame 1`. At the top right, we see that the proportion of reads is greater near the `5'` end of the genes than the `3'` end. Bottom left panel shows the distribution of read lengths in our dataset, while the bottom right panel describes the fraction of bases at the `5'` and `3'` terminals of the reads.

<img src=/examples/reading_frame.png width=350 height=250><img src=/examples/asymmetry.png width=350 height=250>

<img src=/examples/read_dist.png width=350 height=250><img src=/examples/basefrac.png width=350 height=250>

#### P-site offsets
Metagene plots showing normalized ribosome densities for a particular read length (`K-mer`). The well-defined peaks reflect the `3'` boundary of initiating ribosomes as they wait with their p-sites positioned on the start codon (position 0). This allows us to obtain the distance of the ribosome's p-site from the `3'` extrimity of a read, which, as we see here, is a constant `15 nt` for all read lengths.


<img src=/examples/meta_rpf.png width=700 height=500>

#### Initiation
Having the p-site offsets allows us to obtain the metagene profile of all initiation events, averaged across the translatome. We can observe a characteristic three-nucleotide periodicity pattern as the ribosomes wait to accept the cognate tRNA at each codon. Also, since initiation is a rate limiting step in the translation process, the ribosomes spend more time at the start codon on average. This is reflected by the sharp signal at position 0:

<img src=/examples/initiation_peak.png width=700 height=200>
<br></br>

We see how the reading frame is being maintained by the translating ribosomes, with a majority of time being spent on frame 1 (red):

<img src=/examples/start_framing.png width=700 height=200>

#### Pause events
Heatmap reflecting the amount of time a ribosome takes to translate the codon of a particular amino acid. For example, here we see that the ribosomes spend a longer time with histidine codons in their A-sites. The situation is reversed for asparagine or phenylalanine codons.


<img src=/examples/aa_pause_hmap.png width=700 height=250>

<br></br>
We can also break down the above heatmap on a per-codon basis, which shows that both the histidine codons CAT and CAC contribute significantly to the A-site pause. However, only the cognate CAT codon induces a similar pause at the P-site.

<img src=/examples/codon_pause_hmap.png width=700 height=250>

<br></br>

### Processing multiple BAM files in parallel
If you wish to analyze multiple BAM files simultaneously, you can do so by placing all the corresponding configuration files in a single folder, then executing the following command:

```shell
~$ for config in *.json; do echo python bard.py $config; done > jobs.txt

```

This will create a `jobs.txt` file with the following contents:

```bash
python bard.py config_1.json # for the first BAM file
python bard.py config_2.json # and so on ..

...

python bard.py config_N.json
```

You can now use GNU Parallel to process them simultaneously, like so:

```shell
~$ parallel < jobs.txt
```

## License
`bard` is licensed under the GPLv3.
