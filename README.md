#  bard
Batch Analyser of Ribosome profiling datasets

## Requirements

* Python 3.x
* Linux / macOS
* Do a _'pip install requirements.txt'_
* Done!

add info on how to process the gtf before use

bard runs successfully on Linux / macOS.

To get started, you will need the following:

* **Read alignment file in BAM format**
   You need to generate a BAM file


* **Annotation file in GFFv3/GTF format**

* **CDS sequences of the organism in FASTA format**

## Quick example (Yeast)

Enter yeast example images etc here ...


## Installation and usage

bard uses Python 3.x. Please ensure that it is installed and available in your path

First, ensure that you have all the required dependencies installed. Run the following:

```shell
~$ pip install requirements.txt
```

As a user, the **bard.py** file is enough. Run it in your shell like so:

```shell
~$ python bard.py config.json
```

This will write the results in a folder with the following structure:

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
